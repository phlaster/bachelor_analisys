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

# Command-line argument parsing
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

    accessions = readlines(input_file) .|> strip |> unique
    if isempty(accessions)
        @error "Input file '$input_file' is empty"
        exit(1)
    end

    data = DataFrame(accession=String[], gene_count=Int[], genome_length=Int[])
    missing_accessions = String[]

    @showprogress for acc in accessions
        acc_dir = joinpath(genomes_dir, acc, "ncbi_dataset", "data", acc)
        if !isdir(acc_dir)
            push!(missing_accessions, acc)
            continue
        end

        fna_files = filter(f -> endswith(f, "_genomic.fna"), readdir(acc_dir))
        if length(fna_files) != 1
            @warn "Expected one _genomic.fna file in $acc_dir, found $(length(fna_files))"
            continue
        end
        fna_file = joinpath(acc_dir, fna_files[1])

        gff_file = joinpath(acc_dir, "genomic.gff")
        if !isfile(gff_file)
            @warn "GFF file missing for $acc"
            continue
        end

        gene_count = count_genes(gff_file)
        genome_length = get_genome_length(fna_file)
        push!(data, (acc, gene_count, genome_length))
    end

    if !isempty(missing_accessions)
        @warn "Missing directories for $(length(missing_accessions)) accessions: $(join(missing_accessions, ", "))"
    end

    if isempty(data)
        @error "No data to plot"
        exit(1)
    end

    hist = plot(histogram(x=data.gene_count, name="Gene Count"),
        Layout(title="Histogram of Gene Counts",
            xaxis_title="Number of Genes",
            yaxis_title="Frequency"))
    PlotlyJS.savefig(hist, "gene_count_histogram.html")

    scatter_plot = plot(scatter(x=data.genome_length, y=data.gene_count, mode="markers", name="Genome Length vs Gene Count"),
        Layout(title="Genome Length vs Gene Count",
            xaxis_title="Genome Length (bp)",
            yaxis_title="Number of Genes"))
    PlotlyJS.savefig(scatter_plot, "genome_length_vs_gene_count.html")

    @info "Plots saved to gene_count_histogram.html and genome_length_vs_gene_count.html"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end