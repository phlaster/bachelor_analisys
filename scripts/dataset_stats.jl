#!/usr/bin/env julia

const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")
const PLOTS_FILE = joinpath(PROJECT_DIR, "src", "plots.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)
include(PLOTS_FILE)

using ArgParse
using DataFrames
using ProgressMeter
using JSON3
using CSV
using StatsBase

function parse_commandline()
    s = ArgParseSettings(prog="Dataset Statistics",
        description="Get a comprehensive dataset statistics")
    @add_arg_table! s begin
        "--input-file", "-i"
            help = "File with list of accessions"
            required = true
        "--genomes-dir", "-d"
            help = "Directory with individual genomes extracted from NCBI archives"
            required = true
        "--input-train", "-t"
            help = "Optional file with list of accessions for training data"
        "--out-directory", "-o"
            help = "Optional Directory name to save the results"
    end
    return parse_args(s)
end

function check_files_exist(genomes_dir::AbstractString, accession::AbstractString)
    acc_dir = joinpath(genomes_dir, accession)
    fna_file = joinpath(acc_dir, "$accession.fna")
    gff_file = joinpath(acc_dir, "$accession.gff")
    json_file = joinpath(acc_dir, "$accession.json")
    return all(isfile, [fna_file, gff_file, json_file])
end

function get_gene_count(json_file::String)
    json_data = JSON3.read(json_file)
    return json_data.annotationInfo.stats.geneCounts.proteinCoding
end

function count_genes(gff_file::String)
    gff_records = open_gff(gff_file)
    gene_count = count(r -> GFF3.featuretype(r) == "gene", gff_records)
    return gene_count
end

function get_genome_length(fasta_file::String)
    total_length = 0
    open(fasta_file) do f
        for line in eachline(f)
            if !startswith(line, '>')
                total_length += length(strip(line))
            end
        end
    end
    return total_length
end

function process_accessions!(df::DataFrame, genomes_dir::String)
    initial_size = nrow(df)
    changed = false

    filter!(row -> check_files_exist(genomes_dir, row.accession), df)
    if initial_size < nrow(df)
        changed = true
    end
    
    if isempty(df)
        @error "No valid accessions found in $genomes_dir"
        exit(1)
    end

    if !("phylum" in names(df))
        transform!(df, :gtdb_taxonomy => (x -> [split(t, ';')[2][4:end] for t in x]) => :phylum)
        changed = true
    else
        @info "Skip: `phylum` column exists"
    end

    if !("gene_count" in names(df))
        gene_count = Int[]
        @showprogress desc="Collecting gene count from json metadata" for acc in df.accession
            json_file = joinpath(genomes_dir, acc, "$acc.json")
            push!(gene_count, get_gene_count(json_file))
        end
        clear_last_lines(1)
        df.gene_count = gene_count
        changed = true
    else
        @info "Skip: `gene_count` column exists"
    end

    return changed
end

function main()
    args = parse_commandline()
    input_file = abspath(args["input-file"])
    genomes_dir = abspath(args["genomes-dir"])
    train_file = isnothing(args["input-train"]) ? nothing : abspath(args["input-train"])
    out_dir = isnothing(args["out-directory"]) ? nothing : abspath(mkpath(args["out-directory"]))


    if !isfile(input_file)
        @error "Input file '$input_file' does not exist"
        exit(1)
    end
    if !isdir(genomes_dir)
        @error "Directory '$genomes_dir' does not exist"
        exit(1)
    end
    if !isnothing(train_file) && !isfile(train_file)
        @error "Training input file '$train_file' does not exist"
        exit(1)
    end

    main_df = CSV.read(input_file, DataFrame)
    changed = process_accessions!(main_df, genomes_dir)

    if changed
        CSV.write(input_file, main_df; delim='\t')
        @info "Overwritten: $input_file"
    else
        @info "No changes to existing $input_file"
    end
    
    train_df = if !isnothing(train_file)
        train_df = CSV.read(train_file, DataFrame)
        train_changed = process_accessions!(train_df, genomes_dir)
        if train_changed
            CSV.write(train_file, train_df; delim='\t')
            @info "Overwritten: $train_file"
        else
            @info "No changes to existing $train_file"
        end
        train_df
    end
    plotfile = isnothing(out_dir) ? nothing : joinpath(out_dir, "stats.html")
    plt = summaryplots(main_df; train_df=train_df)
    !isnothing(plotfile) && save_plotly(plotfile, plt)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end