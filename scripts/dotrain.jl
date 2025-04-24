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
        "--dump_file", "-d"
            required = true
            help = ""
        "--epoch", "-e"
            default = 3
            arg_type = Int64
            range_tester = x->x≥1
        "--device", "-D"
            default = 0
            arg_type = Int64
            range_tester = x->0≤x≤3
        "--train", "-T"
            default = 1000
            arg_type = Int64
            range_tester = x->x≥1
        "--test", "-t"
            default = 100
            arg_type = Int64
            range_tester = x->x≥1
    end

    return parse_args(s)
end


function main()
    args = parse_commandline()

    N_TRAIN = args["train"]
    N_TEST = args["test"]
    ADDED_EPOCHS = args["epoch"]
    GPU_ID = args["device"]
    
    device!(GPU_ID)
    DEV = gpu_device()
    
    DUMP_FILE = args["dump_file"]
    @assert isfile(DUMP_FILE) "No file $DUMP_FILE"
    dump_dir = dirname(DUMP_FILE)
    dump_data = deserialize(DUMP_FILE)

    @show PAD = dump_data.pad
    @show WINDOW = 2PAD + 1
    @show CHUNK_SIZE = dump_data.chunk
    @show DECAY = dump_data.decay
    @show GAMMA = dump_data.gamma
    @show CHUNK_SKIP = hasfield(typeof(dump_data), :chunk_skip) ? dump_data.chunk_skip : 0.0
    @show LAMBDA = hasfield(typeof(dump_data), :lambda) ? dump_data.lambda : 0.0
    @show ELAPSED_EPOCHS = dump_data.epoch
    
    pseudomonadota_accs = readdlm("DATA/subsets/pseudomonadota.tsv", '\t')[:, 1] .|> string
    
    train_accs = pseudomonadota_accs[1:N_TRAIN]
    test_accs = pseudomonadota_accs[N_TRAIN+1:N_TRAIN+N_TEST]

    dirs_train = "DATA/genomes/genomes/" .* train_accs
    dirs_test = "DATA/genomes/genomes/" .* test_accs
    
    ds_train = GenomeDataset(dirs_train, cds_side=:starts, pad=PAD)
    ds_test = GenomeDataset(dirs_test, cds_side=:starts, pad=PAD)
    
    @show ds_train
    @show ds_test
    
    
    model = dump_data.model
    opt = dump_data.opt
    losses = dump_data.loss
    lrs = dump_data.lr
    metrics = dump_data.metrics
    cms = dump_data.cm

    @info "Dotraining model..."
    dotrain_model(model, ds_train, ds_test;
        epochs=ADDED_EPOCHS,
        dev=DEV,
        opt=opt,
        elapsed_epochs=ELAPSED_EPOCHS,
        floss_gamma=GAMMA,
        loss_lambda=LAMBDA,
        decay_factor=DECAY,
        max_chunk_size=CHUNK_SIZE,
        chunk_skip_coeff=CHUNK_SKIP,
        losses=losses,
        lrs=lrs,
        metrics=metrics,
        cms=cms,
        savedir=dump_dir
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end