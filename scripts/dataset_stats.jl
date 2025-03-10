#!/usr/bin/env julia

const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)

using ArgParse
using DataFrames
using PlotlyJS
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
        "--append", "-a"
            help = "Append new columns to metadata table"
            action = :store_true
    end
    return parse_args(s)
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

function main()
    args = parse_commandline()
    input_file = abspath(args["input-file"])
    genomes_dir = abspath(args["genomes-dir"])

    if !isfile(input_file)
        @error "Input file '$input_file' does not exist"
        exit(1)
    end
    if !isdir(genomes_dir)
        @error "Directory '$genomes_dir' does not exist"
        exit(1)
    end

    df = CSV.read(input_file, DataFrame)
    transform!(df, :gtdb_taxonomy => (x -> [split(t, ';')[2][4:end] for t in x]) => :philum)
    gene_count = Union{Int, Missing}[]

    for (i, acc) in enumerate(df.accession)
        data_dir = joinpath(genomes_dir, acc, "ncbi_dataset", "data")
        acc_dir = joinpath(data_dir, acc)
        
        if !isdir(acc_dir)
            @error "No $acc_dir found"
            push!(gene_count, missing)
            continue
        end

        fna_files = filter(f -> endswith(f, "_genomic.fna"), readdir(acc_dir))
        if length(fna_files) != 1
            @warn "Expected one _genomic.fna file in $acc_dir, found $(length(fna_files))"
            push!(gene_count, missing)
            continue
        end

        gff_file = joinpath(acc_dir, "genomic.gff")
        if !isfile(gff_file)
            @warn "GFF file missing for $acc"
            push!(gene_count, missing)
            continue
        end

        jsonl_file = joinpath(data_dir, "assembly_data_report.jsonl")
        if !isfile(jsonl_file)
            @warn "Assembly report missing for $acc"
            push!(gene_count, missing)
            continue
        end
        
        n_genes = JSON3.read(jsonl_file).annotationInfo.stats.geneCounts.proteinCoding
        push!(gene_count, n_genes)
    end

    df.gene_count = gene_count
    
    cm_philum = countmap(df.philum)
    n_genomes = nrow(df)
    
    hist = plot(histogram(x=df.gene_count, name="Gene Count"),
        Layout(title="Histogram of Gene Counts: $n_genomes",
        xaxis_title="Number of Genes",
        yaxis_title="Frequency")
    )
    PlotlyJS.savefig(hist, "gene_count_histogram.html")
    
    traces = [
        begin
            sub_df = df[df.philum .== p, :]
            scatter(
                x = sub_df.genome_size,
                y = sub_df.gene_count,
                mode = "markers",
                name = "$(cm_philum[p])($(round(100cm_philum[p]/n_genomes, digits=2))): $p",
                marker = attr(
                    opacity = 0.6
                )
            )
        end for p in getindex.(sort(collect(cm_philum), by=x->x.second, rev=true), 1)
    ]

    layout = Layout(
    title = "Genome Length vs Gene Count: $n_genomes",
    xaxis = attr(
        title = "Genome Length (bp)",
        range = [-1e6, 13e6]
    ),
    yaxis = attr(
        title = "Number of Genes",
        range = [0, 12e3]
    ),
    legend = attr(
        x = 0,
        y = 1
    )
)

    scatter_plot = Plot(traces, layout)
    PlotlyJS.savefig(scatter_plot, "genome_length_vs_gene_count.html")

    @info "Plots saved to gene_count_histogram.html and genome_length_vs_gene_count.html"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end