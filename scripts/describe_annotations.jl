#!/usr/bin/env julia
const PROJECT_DIR = dirname(@__DIR__)
const SRC_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(SRC_FILE)

using ArgParse
using Format
using GFF3
using DataFrames
using PrettyTables
using ProgressMeter
using ColorSchemes

Base.show(io::IO, x::Float64) = print(io, format(x, precision=3))

function colored_table(mtr; header, columns_with_numbers)
    hl = color_cols -> Highlighter(
        (data, i, j) -> begin
             j in color_cols && (data[i, j] isa Number) && (0 <= data[i, j] <= 1)
        end,
        (h, data, i, j) -> begin
             val = data[i, j]
             color = get(ColorSchemes.coolwarm, val, (0, 1))
             return Crayon(foreground = (
                 round(Int, color.r * 255),
                 round(Int, color.g * 255),
                 round(Int, color.b * 255)
             ))
        end
    )

    pretty_table(mtr;
        crop = :none,
        tf=tf_borderless,
        header=header,
        highlighters = hl(columns_with_numbers),
        formatters = ft_printf("%.3f")
    )
end


function parse_commandline()
    s = ArgParseSettings()
    s.usage = "./scripts/describe_annotations.jl"
    s.description = "Describe genome annotations in chosen dir across the whole genome"
    s.version = "0.1"
    s.add_version = true

    @add_arg_table! s begin
        "-r", "--ref"
            help = "Reference .gff file (ground truth annotation)"
            required = true
            arg_type = String
        "-d", "--dir"
            help = "Directory where `.gff` annotation files are stored"
            required = true
            arg_type = String
        "-D", "--detailed"
            help = "Print detailed output for each chromosome and strand"
            action = :store_true
    end
    return parse_args(s)
end

function extract_grouped_ranges(filename::String)
    @info "Extracting ranges from $filename"
    records = open_gff(filename)
    chromosomes = unique([GFF3.seqid(r) for r in records])
    strands = ["+", "-"]
    grouped = Dict{Tuple{String, String}, Vector{UnitRange{Int}}}()
    for chrom in chromosomes
        for strand in strands
            filtered = filter_gff_region(sequence_header=chrom, regiontype="CDS", strand=strand, phase=0)(records)
            ranges = [GFF3.seqstart(r):GFF3.seqend(r) for r in filtered]
            grouped[(chrom, strand)] = ranges
        end
    end
    clear_last_lines(1)
    return grouped
end

function compute_scores(TP, FP, FN)
    prec = (TP + FP) > 0 ? TP / (TP + FP) : 0.0
    rec = (TP + FN) > 0 ? TP / (TP + FN) : 0.0
    f1 = (prec + rec) > 0 ? 2 * prec * rec / (prec + rec) : 0.0
    fdr = (TP + FP) > 0 ? FP / (TP + FP) : 0.0
    return (prec=prec, rec=rec, f1=f1, fdr=fdr)
end

function compute_overall_scores(reference_grouped, helixer_grouped, helixer_files, ref_keys)
    total_CMs = @NamedTuple{TP::Int64, FP::Int64, FN::Int64}[]
    for hg in helixer_grouped
        total_TP = 0
        total_FP = 0
        total_FN = 0
        for key in ref_keys
            ref_ranges = get(reference_grouped, key, UnitRange{Int}[])
            pred_ranges = get(hg, key, UnitRange{Int}[])
            cm = ConfusionMTR(ref_ranges, pred_ranges)
            total_TP += cm.TP
            total_FP += cm.FP
            total_FN += cm.FN
        end
        push!(total_CMs, (TP=total_TP, FP=total_FP, FN=total_FN))
    end

    total_TP_comb = 0
    total_FP_comb = 0
    total_FN_comb = 0
    for key in ref_keys
        ref_ranges = get(reference_grouped, key, UnitRange{Int}[])
        all_pred_ranges = UnitRange{Int}[]
        for hg in helixer_grouped
            append!(all_pred_ranges, get(hg, key, UnitRange{Int}[]))
        end
        merged_ranges = merge_ranges(all_pred_ranges)
        cm = ConfusionMTR(ref_ranges, merged_ranges)
        total_TP_comb += cm.TP
        total_FP_comb += cm.FP
        total_FN_comb += cm.FN
    end

    data = Vector{Union{String, Float64}}[]
    for (i, cm) in enumerate(total_CMs)
        scores = compute_scores(cm.TP, cm.FP, cm.FN)
        sample_name = first(split(basename(helixer_files[i]), '.'))
        push!(data, [string(sample_name), scores.prec, scores.rec, scores.f1, scores.fdr])
    end
    scores_comb = compute_scores(total_TP_comb, total_FP_comb, total_FN_comb)
    push!(data, ["combined", scores_comb.prec, scores_comb.rec, scores_comb.f1, scores_comb.fdr])

    return data
end

function add_sample_rows(df, sample_name, reference_grouped, grouped, ref_keys)
    total_TP = 0
    total_FP = 0
    total_FN = 0
    for key in ref_keys
        ref_ranges = get(reference_grouped, key, UnitRange{Int}[])
        pred_ranges = get(grouped, key, UnitRange{Int}[])
        cm = ConfusionMTR(ref_ranges, pred_ranges)
        total_TP += cm.TP
        total_FP += cm.FP
        total_FN += cm.FN
    end
    scores = compute_scores(total_TP, total_FP, total_FN)
    push!(df, (sample_name, "overall", "both", scores.prec, scores.rec, scores.f1, scores.fdr))

    for (chrom, strand) in ref_keys
        ref_ranges = get(reference_grouped, (chrom, strand), UnitRange{Int}[])
        pred_ranges = get(grouped, (chrom, strand), UnitRange{Int}[])
        cm = ConfusionMTR(ref_ranges, pred_ranges)
        scores = compute_scores(cm.TP, cm.FP, cm.FN)
        push!(df, (sample_name, chrom, strand, scores.prec, scores.rec, scores.f1, scores.fdr))
    end
end

function main()
    args = parse_commandline()
    reference_file = args["ref"]
    helixer_dir = args["dir"]
    detailed = args["detailed"]

    helixer_files = filter(endswith(r".[gG][fF][fF]3?"), readdir(helixer_dir, join=true))
    if isempty(helixer_files)
        println("No annotation files found in $helixer_dir\nOnly `.gff` files are accepted.")
        exit()
    end

    reference_grouped = extract_grouped_ranges(reference_file)
    helixer_grouped = [extract_grouped_ranges(file) for file in helixer_files]
    ref_keys = collect(keys(reference_grouped))
    sort!(ref_keys)  # Ensure consistent order

    if detailed
        df = DataFrame(
            Sample=String[],
            Chromosome=String[],
            Strand=String[],
            Precision=Float64[],
            Recall=Float64[],
            F1=Float64[],
            FDR=Float64[]
        )

        @showprogress for (i, file) in enumerate(helixer_files)
            sample_name = first(split(basename(file), '.'))
            add_sample_rows(df, sample_name, reference_grouped, helixer_grouped[i], ref_keys)
        end
        clear_last_lines(1)

        combined_grouped = Dict()
        for key in ref_keys
            all_pred_ranges = UnitRange{Int}[]
            for hg in helixer_grouped
                append!(all_pred_ranges, get(hg, key, UnitRange{Int}[]))
            end
            merged_ranges = merge_ranges(all_pred_ranges)
            combined_grouped[key] = merged_ranges
        end
        add_sample_rows(df, "combined", reference_grouped, combined_grouped, ref_keys)

        println("Detailed annotation scores:")
        colored_table(Matrix(df);
            header=names(df),
            columns_with_numbers=(4,5,6,7)
        )
    else
        data = compute_overall_scores(reference_grouped, helixer_grouped, helixer_files, ref_keys)
        println("Overall annotation scores:")
        colored_table(permutedims(hcat(data...));
            header=["Sample", "Precision", "Recall", "F1", "FDR"],
            columns_with_numbers=(2,3,4,5)
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end