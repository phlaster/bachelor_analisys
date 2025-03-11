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
using CSV
using Downloads
using GZip
using Printf
using JSON3
using StatsBase
using Chain

const METADATA_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz"
const TREE_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_r220.tree.gz"
const GENBANK_PATTERN = r"GC[AF]_\d+\.\d+"
extract_accession(acc) = match(GENBANK_PATTERN, acc).match


function parse_commandline()
    s = ArgParseSettings(
        prog="cluster_genomes.jl",
        description = "Script for preparing dataset accessions through filtering and tree clustering",
        usage="""
        cluster_genomes.jl -i INPUT-DIR [-o OUTPUT-DIR] [-t TRESHOLD] [-m CLUSTER-METHOD] [-p COMPLETENESS-from [- CONTAMINATION-MAX] [-g CONTIGS-MAX] [-r TRNA-MIN] [-n N50-MIN] [-h]

        1. Fetch GTDB metadata+tree via hardcoded links into --input-dir (or read them from there, if they are there);
        2. Cluster the tree using TreeCluster.py with chosen algorithm and threshold;
        3. Filter out accessions, that are not present in the tree file
        4. Filter out using c.l.arg. filters ∪ ncbi_rrna_count≠"none" ∪ ncbi_assembly_level∈{"Complete Genome", "Chromosome"}
        5. Rename accessions to match GenBank naming (eg.: RS_GCF_000657795.2 -> GCF_000657795.2)
        6. Resulting files in `--output-dir`:
            - `cg_sankey.html`: PlotlyJS page visualising filtering steps as https://en.wikipedia.org/wiki/Sankey_diagram
            - `cg_report.json`: report with used parameters, intermediate sample sizes and final stats
            - `cg_hq.tsv`: table with remaining high quality accessions and their relevant fields
            - `cg_train.tsv`: subset of filtered data, where from each cluster a single accession is taken
            - `cg_valid.tsv`: subset of filtered {data\\train}, where from each cluster a single accession is taken
        """)

    @add_arg_table! s begin
        "--input-dir", "-i"
            help = "Output directory for tree and metadata files (required), if no files there, download them using hardcoded links"
            required = true
        "--output-dir", "-o"
            help = "Output file name for script output. --input-dir if not provided"

        "--treshold", "-t"
            help = "Clustering threshold for TreeCluster.py"
            default = 0.2463
            arg_type = Float64
        "--cluster-method", "-m"
            help = "Clustering method (see: https://github.com/niemasd/TreeCluster)"
            default = "avg_clade"
            range_tester = x->in(x, [
                "avg_clade", "leaf_dist_avg", "leaf_dist_max",
                "leaf_dist_min", "length", "length_clade",
                "max", "max_clade", "med_clade", "root_dist",
                "single_linkage", "single_linkage_cut",
                "single_linkage_union", "sum_branch", "sum_branch_clade"])

        "--completeness-min", "-p"
            help = "Minimum completeness percentage (from metadata)"
            default = 98.0
            arg_type = Float64
            range_tester = x->0.0≤x≤100.0
        "--contamination-max", "-c"
            help = "Maximum contamination (from metadata)"
            default = 5.0
            arg_type = Float64
            range_tester = x->0.0≤x≤100.0
        "--contigs-max", "-g"
            help = "Maximum contig count for a genome (from metadata)"
            default = 200
            arg_type = Int64
            range_tester = x->x≥1
        "--trna-min", "-r"
            help = "Minimum tRNA count for a genome (from metadata)"
            default = 40
            arg_type = Int64
            range_tester = x->x≥0
        "--n50-min", "-n"
            help = "N50 score for a genome (from metadata)"
            default = 50_000
            arg_type = Int64
            range_tester = x->x≥1
    end

    return parse_args(s)
end


function download_file(url, dest)
    if !isfile(dest)
        @info "Downloading $url..."
        gz_file = tempname()
        Downloads.download(url, gz_file)
        open(dest, "w") do out
            GZip.open(gz_file) do inp
                write(out, read(inp))
            end
        end
        rm(gz_file)
    end
end

function pick_one_per_cluster(df::DataFrame; skip=String[])
    df_filtered = filter(row -> !(row.accession in skip), df)
    groups = groupby(df_filtered, :cluster)
    picks = [first(g.accession) for g in groups]
    return picks
end


function main()
    args = parse_commandline()
    check_dependencies(["TreeCluster.py", "nw_prune"])

    dir_in = abspath(args["input-dir"])
    dir_out = isnothing(args["output-dir"]) ? dir_in : args["output-dir"]

    threshold = args["treshold"]
    method = args["cluster-method"]

    completeness_min = args["completeness-min"]
    contamination_max = args["contamination-max"]
    contig_count_max = args["contigs-max"]
    trna_count_min = args["trna-min"]
    n50_min = args["n50-min"]
    
    mkpath(dir_in)
    mkpath(dir_out)

    # all output filenames
    output_sankey_filter = joinpath(dir_out, "cg_sankey.html")
    output_json = joinpath(dir_out, "cg_report.json")
    output_all_accessions = joinpath(dir_out, "cg_hq.tsv")
    train_filename = joinpath(dir_out, "cg_train.tsv")
    validation_filename = joinpath(dir_out, "cg_valid.tsv")
    
    metadata_file = joinpath(dir_in, "bac120_metadata_r220.tsv")
    tree_file = joinpath(dir_in, "bac120_r220.tree")

    download_file(METADATA_URL, metadata_file)
    download_file(TREE_URL, tree_file)
    
    cluster_file = tempname()
    @info """Clustering GTDB tree:
        method = $method
        thresh = $threshold
    """
    run(`TreeCluster.py -i $tree_file -o $cluster_file -m $method -t $threshold`)
    cluster_df = DataFrame(CSV.File(cluster_file, delim='\t', header=true), ["accession", "cluster"])
    transform!(cluster_df, :accession => ByRow(extract_accession) => :accession)
    GTDB_tree_leaves = nrow(cluster_df)
    @info "GTDB tree leaves count: $GTDB_tree_leaves"


    metadata_df = DataFrame(CSV.File(metadata_file, delim='\t', header=true, ignorerepeated=true))
    gtdb_genomes_count = nrow(metadata_df)
    @info "GTDB genomes count: $gtdb_genomes_count"

    @info """
    Filtering genome metadata:
        completeness ≥ $completeness_min%
        contamination ≤ $contamination_max%
        contig count ≤ $contig_count_max
        tRNA count ≥ $trna_count_min
        N50 score ≥ $n50_min
    """
    filter_steps = Pair{Int64, String}[]
    joined_df = @chain metadata_df begin
        @aside push!(filter_steps, nrow(_)=>"In metadata file")

        transform!(_, :accession => ByRow(extract_accession) => :accession)
    
        innerjoin(_, cluster_df, on=:accession)
        @aside push!(filter_steps, nrow(_)=>"In tree file")
        
        filter(r->r.ncbi_assembly_level in ["Complete Genome", "Chromosome"], _)
        @aside push!(filter_steps, nrow(_)=>"Genome+Chromosome")
        
        filter(r->r.ncbi_rrna_count != "none", _)
        @aside push!(filter_steps, nrow(_)=>">0 rRNA")
        
        filter(r->min(r.checkm2_completeness, r.checkm_completeness) ≥ completeness_min, _)
        @aside push!(filter_steps, nrow(_)=>"completeness>$completeness_min")
    
        filter(r->max(r.checkm2_contamination, r.checkm_contamination) ≤ contamination_max, _)
        @aside push!(filter_steps, nrow(_)=>"contamination<$contamination_max")
        
        filter(r->r.contig_count ≤ contig_count_max, _)
        @aside push!(filter_steps, nrow(_)=>"contigs≤$contig_count_max")
        
        filter(r->r.trna_count ≥ trna_count_min, _)
        @aside push!(filter_steps, nrow(_)=>"tRNA≥$trna_count_min")
        
        filter(r->r.n50_contigs ≥ n50_min, _)
        @aside push!(filter_steps, nrow(_)=>"n50≥$n50_min")
    end
    create_sankey(filter_steps; savefile=output_sankey_filter)
    
    hq_count = nrow(joined_df)
    @info "Genome count after all filters: $hq_count"

    clusters_cm = countmap(joined_df.cluster)
    n_clusters = haskey(clusters_cm, -1) ? (length(clusters_cm) - 1 + clusters_cm[-1]) : length(clusters_cm)
    singletons = count(==(1), clusters_cm.vals) + get(clusters_cm, -1, 0)
    doubletons = count(==(2), clusters_cm.vals)
    
    @info "Clusters count: $n_clusters"
    chao1 = n_clusters + singletons*(singletons-1)/(2doubletons+1) |> round |> Int
    @info "Chao1 = $chao1"
    
    filtered_hq_accessions = Dict(
        "filter_parameters" => Dict(
            "completeness_min" => completeness_min,
            "contamination_max" => contamination_max,
            "contig_count_max" => contig_count_max,
            "trna_count_min" => trna_count_min,
            "n50_min" => n50_min,
            "TreeCluster_method" => method,
            "TreeCluster_threshold" => threshold,
        ),
        "descriptive" => Dict(
            "gtdb_genomes" => gtdb_genomes_count,
            "genomes_after_filtering" => hq_count,
            "gtdb_tree_leaves" => GTDB_tree_leaves,
            "n_clusters" => n_clusters,
            "singletons" => singletons,
            "doubletons" => doubletons,
            "chao1" => chao1,            
        ),
        "sankey_stages" => collect(enumerate(filter_steps))
    )

    open(output_json, "w") do io
        JSON3.pretty(io, filtered_hq_accessions, JSON3.AlignmentContext(alignment=:Left, indent=2))
    end
    @info "Filtering report saved to: $output_json"

    selected_cols = [
        :accession,
        :genome_size,
        :trna_count,
        :contig_count,
        :n50_contigs,
        :ncbi_rrna_count,
        :gc_percentage,
        :cluster,
        :gtdb_taxonomy,
        :ncbi_assembly_level,
    ]

    
    train_accs = pick_one_per_cluster(joined_df)
    validation_accs = pick_one_per_cluster(joined_df; skip=train_accs)
    
    train_df = filter(row -> row.accession in train_accs, joined_df)
    validation_df = filter(row -> row.accession in validation_accs, joined_df)
    
    CSV.write(output_all_accessions, joined_df[!, selected_cols]; delim='\t')
    @info "Training acccessions saved to $train_filename"
    CSV.write(train_filename, train_df[!, selected_cols]; delim='\t')
    @info "Validation acccessions saved to $validation_filename"
    CSV.write(validation_filename, validation_df[!, selected_cols]; delim='\t')
    @info "All filtered accession saved to $output_all_accessions"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end