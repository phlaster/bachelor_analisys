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


using Flux, Statistics, ProgressMeter
using Flux: onehotbatch, logitcrossentropy, DataLoader, onecold

function create_model(; window_size::Int, input_channels::Int=5, n_classes::Int=4)
    Chain(
        Conv((window_size,), input_channels => 16, pad=0, relu),
        Dropout(0.2),
        Conv((3,), 16 => 32, pad=1, relu),
        Dropout(0.2),
        Conv((1,), 32 => n_classes)
    )
end

function weighted_logitcrossentropy(ŷ::CuArray{Float32, 2}, y::CuArray{Bool, 2}, cw::CuArray{Float32, 1})
    ls = Flux.logsoftmax(ŷ, dims=2)
    W = repeat(reshape(cw, 1, :), size(ŷ, 1), 1)
    loss = -sum(W .* ls .* y) / size(ŷ, 1)
    return loss
end

function train_model!(model, dataset::GenomeDataset; epochs=3, lr=0.001, device=gpu)
    model = model |> device
    opt_state = Flux.setup(Adam(lr), model)
    epoch_losses = Float32[]

    for epoch in 1:epochs
        batch_losses = Float32[]
        @showprogress desc="Epoch $epoch/$epochs" for (i, (seq, labels, weights)) in enumerate(dataset)
            print("$i ")
            seq_dims = size(seq)
            X = reshape(permutedims(seq, (2, 1)), seq_dims[2], seq_dims[1], 1) |> device
            y = permutedims(labels, (2, 1)) |> device
            class_weights = weights |> device

            loss, grads = Flux.withgradient(model) do m
                ŷ = dropdims(m(X), dims=3)
                weighted_logitcrossentropy(ŷ, y, class_weights)
            end
            Flux.update!(opt_state, model, grads[1])
            push!(batch_losses, loss)

        end
        epoch_loss = mean(batch_losses)
        push!(epoch_losses, epoch_loss)
    end
    return epoch_losses
end

function evaluate_model(model, dataset; device=gpu)
    model = model |> device
    n_classes = size(first(dataset)[2], 1)  # Determine number of classes from labels
    true_positives = zeros(Int, n_classes)
    total_per_class = zeros(Int, n_classes)
    
    @showprogress desc="Evaluating" for (seq, labels, _) in dataset
        seq_dims = size(seq)
        X = reshape(permutedims(seq, (2, 1)), seq_dims[2], seq_dims[1], 1) |> device
        y = permutedims(labels, (2, 1)) |> device
        
        ŷ = dropdims(model(X), dims=3)
        preds = onecold(ŷ)
        y_true = onecold(y)
        
        for c in 1:n_classes
            tp = sum((preds .== c) .& (y_true .== c))
            total = sum(y_true .== c)
            true_positives[c] += tp
            total_per_class[c] += total
        end
    end
    
    per_class_accuracy = zeros(n_classes)
    valid_classes = 0
    sum_acc = 0.0
    
    for c in 1:n_classes
        if total_per_class[c] == 0
            per_class_accuracy[c] = NaN
        else
            acc = true_positives[c] / total_per_class[c]
            per_class_accuracy[c] = acc
            sum_acc += acc
            valid_classes += 1
        end
    end
    
    mean_accuracy = valid_classes == 0 ? NaN : sum_acc / valid_classes
    return per_class_accuracy, mean_accuracy
end


function main()
    dirs = "DATA/genomes/genomes/" .* readdir("DATA/genomes/genomes")
    train_dirs, val_dirs = dirs[1:10], dirs[1:10]

    PAD = 15
    WINDOW = 2PAD + 1
    
    ds_train = GenomeDataset(train_dirs, cds_side=:starts, pad=PAD)
    ds_val = GenomeDataset(val_dirs, cds_side=:starts, pad=PAD)

    device = gpu_device()

    model = create_model(; window_size=WINDOW)

    losses = train_model!(model, ds_train; epochs=3, lr=0.01, device=device)
    evaluate_model(model, ds_val; device=device)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end