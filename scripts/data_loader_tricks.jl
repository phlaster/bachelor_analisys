#!/usr/bin/env julia

const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)

using LuxCUDA
using Flux
using MLUtils
using Statistics
using ProgressMeter
using Functors
using Serialization
using DelimitedFiles

function _remove_last_dim(x)
    dropdims(x, dims=3)
end

function create_model(; window_size::Int, input_channels::Int=5, n_classes::Int=4)
    Chain(
        Conv((window_size,), input_channels => 32, pad=0, bias=false),
        BatchNorm(32, relu),
        Dropout(0.3),

        Conv((3,), 32 => 64, pad=1, bias=false),
        BatchNorm(64, relu),
        Dropout(0.3),
        Conv((3,), 64 => 64, pad=1, bias=false),
        BatchNorm(64, relu),
        Dropout(0.3),

        Conv((1,), 64 => n_classes),
        _remove_last_dim
    )
end

function weighted_logitcrossentropy(ŷ::CuArray{Float32,2}, y::CuArray{Bool,2}, cw::CuArray{Float32,1})
    ls = Flux.logsoftmax(ŷ, dims=2)
    W = repeat(reshape(cw, 1, :), size(ŷ, 1), 1)
    loss = -sum(W .* ls .* y) / size(ŷ, 1)
    return loss
end

_transform_X(seq) = reshape(permutedims(seq, (2, 1)), size(seq, 2), size(seq, 1), 1)
_transform_y(labels) = permutedims(labels, (2, 1))

function split_into_chunks(seq, labels, max_chunk_size)
    seq_length = size(seq, 2)
    labels_length = size(labels, 2)
    pad = (seq_length - labels_length) ÷ 2
    seq_length <= max_chunk_size && return [(seq, labels)]

    nchunks = labels_length ÷ max_chunk_size + 1
    chunk_lengths_approx = ceil(Int, labels_length/nchunks)
    chunk_coords = Base.Iterators.partition(1:labels_length, chunk_lengths_approx)

    res = @views (
        (seq[:, (coords.start):(coords.stop+2pad)], labels[:, coords])
        for coords in chunk_coords
    )
   return res 
end

function train_model!(model, dataset::GenomeDataset; epochs=1:3, lr=0.001, device=gpu, max_chunk_size=5*10^4, save_prefix="")
    function _add_grads(a, b)
        if isnothing(a)
            if isnothing(b)
                return nothing
            else
                return b
            end
        elseif isnothing(b)
            return a
        end
        return a .+ b
    end
    _scale_grads(x, factor) = isnothing(x) ? nothing : x .* factor
    
    model = model |> device
    opt_state = Flux.setup(Adam(lr), model)
    epoch_losses = Float32[]
    epoch_loss = Inf

    for epoch in epochs
        batch_losses = Float32[]
        @showprogress desc="Epoch $epoch/$(epochs.stop)" for (seq, labels) in dataset
            accumulated_grads = nothing
            accum_counter = 0

            for (seq_chunk, labels_chunk) in split_into_chunks(seq, labels, max_chunk_size)
                X = _transform_X(seq_chunk) |> device
                _y = _transform_y(labels_chunk)
                class_weights = dropdims(sum(_y, dims=1) .* inv(size(_y, 1)), dims=1) |> device
                y = _y |> device

                loss, grads = Flux.withgradient(model) do m
                    ŷ = m(X)
                    weighted_logitcrossentropy(ŷ, y, class_weights)
                end

                accumulated_grads = Functors.fmap(_add_grads, accumulated_grads, grads[1])
                accum_counter += 1
                push!(batch_losses, loss)
            end

            if accum_counter > 0
                scaled_grads = Functors.fmap(x -> _scale_grads(x, inv(accum_counter)), accumulated_grads)
                Flux.update!(opt_state, model, scaled_grads)
            end
        end

        mean_batch_loss = mean(batch_losses)
        if mean_batch_loss < epoch_loss
            serialize("DATA/saved_models/$save_prefix.epoch=$epoch.flux", model)
        end
        epoch_loss = mean(batch_losses)
        
        @info "Mean loss for epoch $epoch: $epoch_loss"
        push!(epoch_losses, epoch_loss)
    end
    return epoch_losses
end

function evaluate_model(model, dataset; device=gpu, max_chunk_size=5*10^4)
    model = model |> device
    n_classes = size(first(dataset)[2], 1)

    cm_classes = [zeros(Int, 2,2) for _ in 1:n_classes]

    @showprogress desc = "Evaluating" for (seq, labels) in dataset
        for (seq_chunk, labels_chunk) in split_into_chunks(seq, labels, max_chunk_size)        
            X = _transform_X(seq_chunk)    |> device
            y_pred = softmax(model(X), dims=2) |> cpu
            predicted_classes = argmax.(eachrow(y_pred))
            true_classes = Flux.onecold(labels_chunk, 1:4)
            cm_classes .+= compute_confusion_matrices(true_classes, predicted_classes)
        end
    end

    return cm_classes
end

function main()
    PAD = 500
    WINDOW = 2PAD + 1
    TOP_N = 1000
    n_epochs = 10

    pseudomonadota_accs = readdlm("DATA/subsets/pseudomonadota.tsv", '\t')[:, 1] .|> string
    SAVEPREFIX = "psdomndt_top1000_starts_PAD=$PAD"
    dirs = "DATA/genomes/genomes/" .* pseudomonadota_accs[1:TOP_N]

    ds_train = GenomeDataset(dirs, cds_side=:starts, pad=PAD,  max_cache_entries=TOP_N)
    @show ds_train
    device = gpu_device()
    model = if isempty(ARGS)
        create_model(; window_size=WINDOW)
    else 
        _mo = deserialize(ARGS[1])
        model_fixed = Chain(_mo.layers[1:end-1]..., _remove_last_dim)
    end

    epochs = if !isempty(ARGS)
        e_start = parse(Int, last(split(split(ARGS[1], '.')[end-1], '=')))
        e_start+1:e_start+n_epochs
    else
        1:n_epochs
    end
    @show epochs

    losses = train_model!(model, ds_train; epochs=epochs, lr=0.001, device=device, save_prefix=SAVEPREFIX)
    conf_matrixes = evaluate_model(model, ds_train; device=device)

    allstats = (losses=losses, conf_mtrx=conf_matrixes)
    serialize("DATA/saved_models/stats_$SAVEPREFIX", allstats)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end