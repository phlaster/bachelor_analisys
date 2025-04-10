#!/usr/bin/env julia

const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)

using DelimitedFiles
using Serialization
using ArgParse
using Dates

function parse_commandline()
    s = ArgParseSettings(
        prog="XXX.jl",
        description = "...",
        usage="""
        """)

    @add_arg_table! s begin
        "--output_dir", "-o"
            required = true
            help = ""
        "--pad", "-p"
            help = ""
            default = 30
            arg_type = Int64
            range_tester = x->x≥0
        "--train", "-T"
            default = 1000
            arg_type = Int64
            range_tester = x->x≥1
        "--test", "-t"
            default = 100
            arg_type = Int64
            range_tester = x->x≥1
        "--epoch", "-e"
            default = 10
            arg_type = Int64
            range_tester = x->x≥1
        "--lr", "-l"
            default = 5e-3
            arg_type = Float64
            range_tester = x->0<x<1
        "--decay", "-d"
            default = 0.6
            arg_type = Float64
            range_tester = x->0<x≤1
        "--gamma", "-g"
            default = 3.0
            arg_type = Float64
            range_tester = x->0≤x
        "--device", "-D"
            default = 0
            arg_type = Int64
            range_tester = x->0≤x≤3
    end

    return parse_args(s)
end


function main()
    args = parse_commandline()
    PAD = args["pad"]
    DIR = args["output_dir"]
    WINDOW = 2PAD + 1
    N_TRAIN = args["train"]
    N_TEST = args["test"]
    n_epochs = args["epoch"]
    lr = args["lr"]
    decay_factor = args["decay"]
    floss_gamma = args["gamma"]
    location = mkpath(joinpath(DIR, string(now())))
    N_GPU = args["device"]
    
    device!(N_GPU)
    

    pseudomonadota_accs = readdlm("DATA/subsets/pseudomonadota.tsv", '\t')[:, 1] .|> string
    
    train_accs = pseudomonadota_accs[1:N_TRAIN]
    test_accs = pseudomonadota_accs[N_TRAIN+1:N_TRAIN+N_TEST]

    dirs_train = "DATA/genomes/genomes/" .* train_accs
    dirs_test = "DATA/genomes/genomes/" .* test_accs
    
    ds_train = GenomeDataset(dirs_train, cds_side=:starts, pad=PAD,  max_cache_entries=N_TRAIN)
    ds_test = GenomeDataset(dirs_test, cds_side=:starts, pad=PAD,  max_cache_entries=N_TEST)
    
    @show ds_train
    @show ds_test
    
    dev = gpu_device()
    
    model = create_model(; window_size=WINDOW)

    model, losses, lrs = train_model(model, ds_train, ds_test;
        epochs=n_epochs,
        lr=lr,
        dev=dev,
        floss_gamma=floss_gamma,
        decay_factor=decay_factor,
        savedir=location
    )

end

if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end