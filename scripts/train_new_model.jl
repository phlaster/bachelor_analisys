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
            default = "."
            help = ""
        "--suffix", "-s"
            default = ""
            help = ""
        "--pad", "-p"
            help = ""
            default = 125
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
            default = 2e-3
            arg_type = Float64
            range_tester = x->0<x<1
        "--decay", "-d"
            default = 0.8
            arg_type = Float64
            range_tester = x->0<x≤1
        "--gamma", "-g"
            default = 3.5
            arg_type = Float64
            range_tester = x->0≤x
        "--lambda", "-b"
            default = 1.f-6
            arg_type = Float32
            range_tester = x->0≤x
        "--skip_chunks", "-x"
            default = 0.98
            arg_type = Float64
            range_tester = x->x>0
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
    WINDOW = 2PAD + 1
    N_TRAIN = args["train"]
    N_TEST = args["test"]
    n_epochs = args["epoch"]
    lr = args["lr"]
    decay_factor = args["decay"]
    floss_gamma = args["gamma"]
    loss_lambda = args["lambda"]
    chunk_skip_coeff = args["skip_chunks"]
    
    DIR = args["output_dir"]
    DIR_SUFFIX = args["suffix"]
    thistime = string(now())
    location = joinpath(DIR, thistime*DIR_SUFFIX)
    N_GPU = args["device"]
    
    for (k, v) in args
        println(rpad(k, 15, " "), ": ", v)
    end
    println(rpad("timestamp", 15, " "), ": ", thistime)

    device!(N_GPU)
    

    pseudomonadota_accs = readdlm("DATA/subsets/pseudomonadota.tsv", '\t')[:, 1] .|> string
    
    train_accs = pseudomonadota_accs[1:N_TRAIN]
    test_accs = pseudomonadota_accs[N_TRAIN+1:N_TRAIN+N_TEST]

    dirs_train = "DATA/genomes/genomes/" .* train_accs
    dirs_test = "DATA/genomes/genomes/" .* test_accs
    
    ds_train = GenomeDataset(dirs_train, cds_side=:starts, pad=PAD)
    ds_test = GenomeDataset(dirs_test, cds_side=:starts, pad=PAD)
    
    @show ds_train
    @show ds_test
    
    dev = gpu_device()
    
    model = create_model(; window_size=WINDOW)

    train_model(model, ds_train, ds_test;
        epochs=n_epochs,
        lr=lr,
        dev=dev,
        floss_gamma=floss_gamma,
        loss_lambda=loss_lambda,
        decay_factor=decay_factor,
        chunk_skip_coeff=chunk_skip_coeff,
        savedir=location
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end