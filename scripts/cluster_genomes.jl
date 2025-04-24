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
using ProgressMeter

const METADATA_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz"
const TREE_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_r220.tree.gz"
const CLUSTERS_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/auxillary_files/sp_clusters_r220.tsv"
const GENBANK_PATTERN = r"GC[AF]_\d+\.\d+"
extract_accession(acc) = match(GENBANK_PATTERN, acc).match


function parse_commandline()
    s = ArgParseSettings(
        prog="cluster_genomes.jl",
        description = "Script for preparing dataset accessions through filtering and tree clustering",
        usage="""
        cluster_genomes.jl -i INPUT-DIR [-o OUTPUT-DIR] [-t TRESHOLD] [-m CLUSTER-METHOD] [-p COMPLETENESS-from [- CONTAMINATION-MAX] [-g CONTIGS-MAX] [-r TRNA-MIN] [-n N50-MIN] [-h]

        1. Fetch GTDB metadata+tree+sp_clusters via hardcoded links into --input-dir (or read them from there, if they are there);
        2. For each leaf in tree file find a representative in same cluster, which has either "Complete Genome" or "Chromosome" status, then write the new substituted tree file.
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
        if endswith(url, ".gz")
            open(dest, "w") do out
                GZip.open(gz_file) do inp
                    write(out, read(inp))
                end
            end
            rm(gz_file)
        else
            mv(gz_file, dest)
        end
    end
end

function parse_tree_leaves(tree_file::String)
    tree_str = read(tree_file, String)
    leaves = Set{String}(m.match for m in eachmatch(GENBANK_PATTERN, tree_str))
    return collect(leaves)
end

function substitute_tree_representatives(tree_file::String, clusters_file::String, metadata_file::String, output_tree_file::String)
    metadata_df = DataFrame(CSV.File(metadata_file, delim='\t', header=true, ignorerepeated=true))
    transform!(metadata_df, :accession => ByRow(extract_accession) => :accession)
    
    clusters_df = DataFrame(CSV.File(clusters_file, delim='\t', header=true))
    transform!(clusters_df, :("Representative genome") => ByRow(extract_accession) => :representative)
    transform!(clusters_df, :("Clustered genomes") => ByRow(x -> split(x, ",")) => :clustered)

    hq_genomes = filter(r -> r.ncbi_assembly_level in ["Complete Genome", "Chromosome"], metadata_df)
    hq_dict = Dict(r.accession => r.ncbi_assembly_level for r in eachrow(hq_genomes))

    sub_dict = Dict{String, String}()
    @showprogress desc = "Build substitution map" for row in eachrow(clusters_df)
        rep = row.representative
        if haskey(hq_dict, rep)
            sub_dict[rep] = rep  # Representative is already high-quality
        else
            # Find a high-quality genome in the cluster
            for acc in row.clustered
                acc_clean = extract_accession(acc)
                if haskey(hq_dict, acc_clean)
                    sub_dict[rep] = acc_clean
                    break
                end
            end
        end
    end

    tree_str = read(tree_file, String)
    for (old, new) in sub_dict
        if old != new
            tree_str = replace(tree_str, old => new)
        end
    end
    write(output_tree_file, tree_str)
    @info "Substituted tree written to $output_tree_file"
    return sub_dict
end

function pick_one_per_cluster(df::DataFrame; skip=String[])
    df_filtered = filter(row -> !(row.accession in skip), df)
    groups = groupby(df_filtered, :cluster)
    picks = String[]
    for g in groups
        if first(g.cluster) == -1
            append!(picks, g.accession)
        else
            push!(picks, first(g.accession))
        end
    end
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

    output_sankey_filter = joinpath(dir_out, "cg_sankey.html")
    output_json = joinpath(dir_out, "cg_report.json")
    output_all_accessions = joinpath(dir_out, "cg_hq.tsv")
    train_filename = joinpath(dir_out, "cg_train.tsv")
    validation_filename = joinpath(dir_out, "cg_valid.tsv")
    
    metadata_file = joinpath(dir_in, "bac120_metadata_r220.tsv")
    tree_file = joinpath(dir_in, "bac120_r220.tree")
    clusters_file = joinpath(dir_in, "sp_clusters_r220.tsv")
    substituted_tree_file = joinpath(dir_in, "bac120_r220_substituted.tree")

    download_file(METADATA_URL, metadata_file)
    download_file(TREE_URL, tree_file)
    download_file(CLUSTERS_URL, clusters_file)

    sub_dict = substitute_tree_representatives(tree_file, clusters_file, metadata_file, substituted_tree_file)
    @info "Substituted $(length(sub_dict)) representatives"

    cluster_file = tempname()
    @info """Clustering substituted GTDB tree:
        method = $method
        thresh = $threshold
    """
    run(`TreeCluster.py -i $substituted_tree_file -o $cluster_file -m $method -t $threshold`)
    cluster_df = DataFrame(CSV.File(cluster_file, delim='\t', header=true), ["accession", "cluster"])
    transform!(cluster_df, :accession => ByRow(extract_accession) => :accession)
    GTDB_tree_leaves = nrow(cluster_df)
    @info "Substituted tree leaves count: $GTDB_tree_leaves"

    metadata_df = DataFrame(CSV.File(metadata_file, delim='\t', header=true, ignorerepeated=true))
    gtdb_genomes_count = nrow(metadata_df)
    @info "GTDB genomes count: $gtdb_genomes_count"

    filter_steps = Pair{Int64, String}[]
    joined_df = @chain metadata_df begin
        @aside push!(filter_steps, nrow(_) => "In metadata file")
        transform!(_, :accession => ByRow(extract_accession) => :accession)
        innerjoin(_, cluster_df, on=:accession)
        @aside push!(filter_steps, nrow(_) => "In substituted tree")
        filter(r -> r.ncbi_assembly_level in ["Complete Genome", "Chromosome"], _)
        @aside push!(filter_steps, nrow(_) => "Genome+Chromosome")
        filter(r -> !(r.ncbi_rrna_count in ["none", "0"]), _)
        @aside push!(filter_steps, nrow(_) => ">0 rRNA")
        filter(r -> min(r.checkm2_completeness, r.checkm_completeness) ≥ completeness_min, _)
        @aside push!(filter_steps, nrow(_) => "completeness>$completeness_min")
        filter(r -> max(r.checkm2_contamination, r.checkm_contamination) ≤ contamination_max, _)
        @aside push!(filter_steps, nrow(_) => "contamination<$contamination_max")
        filter(r -> r.contig_count ≤ contig_count_max, _)
        @aside push!(filter_steps, nrow(_) => "contigs≤$contig_count_max")
        filter(r -> r.trna_count ≥ trna_count_min, _)
        @aside push!(filter_steps, nrow(_) => "tRNA≥$trna_count_min")
        filter(r -> r.n50_contigs ≥ n50_min, _)
        @aside push!(filter_steps, nrow(_) => "n50≥$n50_min")
    end
    sankey_plot = create_sankey(filter_steps)
    save_plotly(output_sankey_filter, sankey_plot)
    
    hq_count = nrow(joined_df)
    @info "Genome count after all filters: $hq_count"

    clusters_cm = countmap(joined_df.cluster)
    n_clusters = haskey(clusters_cm, -1) ? (length(clusters_cm) - 1 + clusters_cm[-1]) : length(clusters_cm)
    singletons = count(==(1), clusters_cm.vals) + get(clusters_cm, -1, 0)
    doubletons = count(==(2), clusters_cm.vals)
    chao1 = n_clusters + singletons * (singletons - 1) / (2 * doubletons + 1) |> round |> Int

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
    CSV.write(train_filename, train_df[!, selected_cols]; delim='\t')
    CSV.write(validation_filename, validation_df[!, selected_cols]; delim='\t')
    @info "Results saved to $output_all_accessions, $train_filename, $validation_filename"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end