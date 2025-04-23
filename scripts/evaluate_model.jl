#!/usr/bin/env julia

const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")
const PLOTS_FILE = joinpath(PROJECT_DIR, "src", "plots.jl")


using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)
include(PLOTS_FILE)


using DelimitedFiles
using Serialization
using ArgParse
using Dates
using Random

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
        "--device", "-D"
            default = 0
            arg_type = Int64
            range_tester = x->0≤x≤3
        "--test", "-t"
            default = 100
            arg_type = Int64
            range_tester = x->x≥1
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()

    N_TEST = args["test"]
    GPU_ID = args["device"]
    
    device!(GPU_ID)
    DEV = gpu_device()
    
    DUMP_FILE = args["dump_file"]
    @assert isfile(DUMP_FILE) "No file $DUMP_FILE"
    DUMP_DIR = dirname(DUMP_FILE)
    dump_data = deserialize(DUMP_FILE)

    @show PAD = dump_data.pad
    @show WINDOW = 2PAD + 1
    
    pseudomonadota_accs = readdlm("DATA/subsets/pseudomonadota.tsv", '\t')[:, 1] .|> string |> shuffle
    test_accs = pseudomonadota_accs[1:N_TEST]
    dirs_test = "DATA/genomes/genomes/" .* test_accs
    ds_test = GenomeDataset(dirs_test, cds_side=:starts, pad=PAD)
    
    @show ds_test
    
    
    model = dump_data.model |> DEV

    (class_metrics, conf_mtr, classes), fp_shifts, gt_s2s_ranges, pred_s2s_ranges = evaluate_bin_class_model(
        model, ds_test; dev=DEV
    )

    serialize(DUMP_FILE * ".stats", (fp_shifts=fp_shifts, gt_s2s_ranges=gt_s2s_ranges, pred_s2s_ranges=pred_s2s_ranges))
    # plot_path = DUMP_FILE * ".png"
    # plot_two_subplots(fp_shifts, gt_s2s_ranges, pred_s2s_ranges; savename=plot_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    @time main()
end