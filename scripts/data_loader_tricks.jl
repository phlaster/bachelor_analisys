#!/usr/bin/env julia

const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)

using DelimitedFiles
using Serialization


function main()
    PAD = 30
    WINDOW = 2PAD + 1
    N_TRAIN = 1000
    N_TEST = 100
    n_epochs = 3
    lr = 5e-3

    pseudomonadota_accs = readdlm("DATA/subsets/pseudomonadota.tsv", '\t')[:, 1] .|> string
    
    dirs_train = "DATA/genomes/genomes/" .* pseudomonadota_accs[1:N_TRAIN]
    ds_train = GenomeDataset(dirs_train, cds_side=:starts, pad=PAD,  max_cache_entries=N_TRAIN)
    @show ds_train
    
    dirs_test = "DATA/genomes/genomes/" .* pseudomonadota_accs[N_TRAIN+1:N_TRAIN+N_TEST]
    ds_test = GenomeDataset(dirs_test, cds_side=:starts, pad=PAD,  max_cache_entries=N_TEST)
    @show ds_test

    device = gpu_device()
    
    model = if isempty(ARGS)
        create_model(; window_size=WINDOW)
    else 
        deserialize(ARGS[1])
    end |> device

    epochs = if !isempty(ARGS)
        e_start = parse(Int, last(split(split(ARGS[1], '.')[end-1], '=')))
        e_start+1:e_start+n_epochs
    else
        1:n_epochs
    end
    
    losses = Float32[]
    lrs = Float64[]
    for _ in 1:4
        @show epochs

        losses_i = train_model!(model, ds_train; epochs=epochs, lr=lr, device=device)
        serialize("DATA/saved_models/BinClasses_model_epochs=$(epochs.stop).flux", model)
        
        cms = evaluate_model(model, ds_test; device=device)
        append!(losses, losses_i)
        push!(lrs, lr)

        allstats = (loss=losses, cm=cms, lr=lrs, pad=PAD, n_train=N_TRAIN, n_test=N_TEST)

        serialize("DATA/saved_models/BinClasses_stats_epochs=$(epochs.stop).bin", allstats)
        
        epochs = (epochs.start+n_epochs):(epochs.stop+n_epochs)
        lr /= 2.5
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end