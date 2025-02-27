#!/usr/bin/env julia
PROJECT_DIR = dirname(@__DIR__)
SRC_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(SRC_FILE)

using ArgParse
using Format
using REPL.TerminalMenus
using GFF3
using DataFrames
using PrettyTables

Base.show(io::IO, x::Float64) = print(io, format( x, precision=3 ))

function parse_commandline()
    s = ArgParseSettings()
    s.usage = "./scripts/describe_annotations.jl"
    s.description = "Describe genome annotations in chosen dir"
    # s.epilog = ""
    s.version = "0.1"
    s.add_version = true

    @add_arg_table! s begin
        "-r", "--ref"
            help = "Filename with reference annotation"
            required = true
            arg_type = String
        "-d", "--dir"
            help = "Directory where `.gff` annotation files are stored"
            required = true
            arg_type = String
        "-m", "--manual"
            help = "Manually choose chromosome for the analisys (alphanum first one by default)"
            action = :store_true
        "-s", "--strand"
            help = "Choose a strand (`+` by default)"
            arg_type = String
            default = "+"
            range_tester = x->in(x, ["+", "-"])
    end
    return parse_args(s)
end


function main()
    args = parse_commandline()
    reference_file = args["ref"]
    helixer_dir = args["dir"]
    
    helixer_files = filter(endswith(".gff"), readdir(helixer_dir, join=true))
    if isempty(helixer_files)
        println("No annotation files found in $helixer_dir\nOnly `.gff` files are accepted.")
        exit()
    end
    chromosomes = Set([GFF3.seqid(x) for x in open_gff(reference_file)]) |> collect |> sort
    
    chrom = if args["manual"]
        menu_height = 4
        menu = RadioMenu(chromosomes, pagesize=menu_height, ctrl_c_interrupt=false)
        menu_choice = request("Choose a chromosome to analyze annotation on:", menu)
        clear_last_lines(menu_height+1)
        menu_choice != -1 ? chromosomes[menu_choice] : (println("No chromosome chosen, exit."); exit())
    else
        first(chromosomes)
    end
    
    pipeline = ranges_from_GFF_records ∘ filter_gff_region(chrom; regiontype="CDS", strand=args["strand"], phase=0) ∘ open_gff
    
    reference_ranges = pipeline(reference_file)
    helixer_ranges = pipeline.(helixer_files)
    ranges_all = merge_ranges(helixer_ranges)

    CMs =  ConfusionMTR.(Ref(reference_ranges), helixer_ranges)
    CM_comb = ConfusionMTR(reference_ranges, ranges_all)

    data = [[
        first(split(basename(name), '.')), cm.prec, cm.rec, cm.f1, cm.fdr
    ] for (name, cm) in zip(helixer_files, CMs)]
	push!(data, ["combined", CM_comb.prec, CM_comb.rec, CM_comb.f1, CM_comb.fdr])
	
    println("Annotation scores for chromosome $chrom:")
	pretty_table(permutedims(hcat(data...)); header=["Sample", "Precision", "Recall", "F1", "FDR"])
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
