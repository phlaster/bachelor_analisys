#!/usr/bin/env julia
const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)

using ArgParse
using DataFrames
using CSV
using Downloads
using GZip
using Printf
using JSON3
using StatsBase

const METADATA_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz"
const TREE_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_r220.tree.gz"
const GENBANK_PATTERN = r"GC[AF]_\d+\.\d+"
extract_accession(acc) = match(GENBANK_PATTERN, acc).match


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--directory", "-d"
            help = "Output directory for tree and metadata files (required)"
            required = true
        "--output", "-o"
            help = "Output file name for selected accessions"

        "--TreeClusterTreshold", "-T"
            help = "Clustering threshold for TreeCluster.py (required)"
            default = 0.6
            arg_type = Float64
        "--TreeClusterMethod", "-m"
            help = "Clustering method"
            default = "avg_clade"
            range_tester = x->in(x, [
                "avg_clade", "leaf_dist_avg", "leaf_dist_max",
                "leaf_dist_min", "length", "length_clade",
                "max", "max_clade", "med_clade", "root_dist",
                "single_linkage", "single_linkage_cut",
                "single_linkage_union", "sum_branch", "sum_branch_clade"])
        
        "--CompletenessMin", "-C"
            help = "Minimum completeness percentage"
            default = 95.0
            arg_type = Float64
            range_tester = x->0.0≤x≤100.0
        "--ContaminationMax", "-c"
            help = "Maximum contamination"
            default = 5.0
            arg_type = Float64
            range_tester = x->0.0≤x≤100.0
        "--ContigCountMax", "-g"
            help = "Maximum contig count for a genome"
            default = 200
            arg_type = Int64
            range_tester = x->x≥1
        "--tRNACountMin", "-R"
            help = "Minimum tRNA count for a genome"
            default = 10
            arg_type = Int64
            range_tester = x->x≥0
        "--N50Min", "-N"
            help = "N50 score for a genome"
            default = 30_000
            arg_type = Int64
            range_tester = x->x≥1
        "--StrainHeterogenityMax", "-H"
            help = "Maximum strain heterogenity score"
            default = 10.0
            arg_type = Float64
            range_tester = x->0.0≤x≤100.0
    end

    return parse_args(s)
end

function check_dependencies()
    try
        run(pipeline(`TreeCluster.py --help`, stdout=devnull, stderr=devnull))
    catch e
        @error "treecluster is not installed or not in PATH"
        exit()
    end
    try
        run(pipeline(`nw_prune -h`, stdout=devnull, stderr=devnull))
    catch e
        @error "newick_utils is not installed or not in PATH"
        exit()
    end
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

function chao1(cm, f₁, f₂)
    S_obs = length(cm)
    return f₂ > 0 ? S_obs + f₁^2/2f₂ : S_obs + f₁*(f₁-1)/2
end


function main()
    check_dependencies()
    args = parse_commandline()

    dir = abspath(args["directory"])
    output = args["output"]

    threshold = args["TreeClusterTreshold"]
    method = args["TreeClusterMethod"]

    completeness_min = args["CompletenessMin"]
    contamination_max = args["ContaminationMax"]
    contig_count_max = args["ContigCountMax"]
    trna_count_min = args["tRNACountMin"]
    n50_min = args["N50Min"]
    sh_max = args["StrainHeterogenityMax"]
    
    
    mkpath(dir)
    
    metadata_file = joinpath(dir, "bac120_metadata_r220.tsv")
    tree_file = joinpath(dir, "bac120_r220.tree")

    download_file(METADATA_URL, metadata_file)
    download_file(TREE_URL, tree_file)

    
    metadata_df = DataFrame(CSV.File(metadata_file, delim='\t', header=true, ignorerepeated=true))
    @info "Metadata file entries: $(nrow(metadata_df))"

    @info """
    Filtering genome metadata:
        completeness ≥ $completeness_min%
        contamination ≤ $contamination_max%
        contig count ≤ $contig_count_max
        tRNA count ≥ $trna_count_min
        N50 score ≥ $n50_min
        Strain Heterogenity ≤ $sh_max%
    """
    hq_metadata_df = filter(row -> 
        min(row.checkm2_completeness, row.checkm_completeness) ≥ completeness_min &&
        max(row.checkm2_contamination, row.checkm_contamination) ≤ contamination_max &&
        row.contig_count ≤ contig_count_max &&
        row.trna_count ≥ trna_count_min &&
        row.n50_contigs ≥ n50_min &&
        row.checkm_strain_heterogeneity ≤ sh_max &&
        
        row.ambiguous_bases == 0 &&
        # occursin.(strip.(first.(split.(row.ncbi_organism_name)), Ref(['[', ']'])), row.gtdb_taxonomy) &&
        row.mimag_high_quality == "t" &&
        row.ncbi_assembly_level in ["Complete Genome", "Chromosome"] &&
        row.ncbi_rrna_count != "none" && true,
        # row.gtdb_type_designation_ncbi_taxa_sources != "none",
    metadata_df)
    transform!(hq_metadata_df, :accession => ByRow(extract_accession) => :accession)
    n_hq = nrow(hq_metadata_df)
    @info "Genome count after filtering: $(n_hq)"
        
    cluster_file = tempname()
    @info """Clustering GTDB tree:
        method=$method
        threshold=$threshold
    """
    run(`TreeCluster.py -i $tree_file -o $cluster_file -m $method -t $threshold`)
    cluster_df = DataFrame(CSV.File(cluster_file, delim='\t', header=true), ["accession", "cluster"])
    transform!(cluster_df, :accession => ByRow(extract_accession) => :accession)
    @info "GTDB tree leaves count: $(nrow(cluster_df))"
    
    joined_df = innerjoin(hq_metadata_df, cluster_df, on=:accession)
    n_clustered_genomes = nrow(joined_df)
    @info "Joined metadata accessions and tree accessions count: $n_clustered_genomes"

    clusters_cm = countmap(joined_df.cluster)
    n_clusters = length(clusters_cm) - 1 + get(clusters_cm, -1, 0)
    singletons = count(==(1), clusters_cm.vals) + get(clusters_cm, -1, 0)
    doubletons = count(==(2), clusters_cm.vals)
    
    @info "Clusters count: $n_clusters"
    
    ###############
    filtered_hq_accessions = Dict(
        "parameters" => Dict(
            "completeness_min" => completeness_min,
            "contamination_max" => contamination_max,
            "contig_count_max" => contig_count_max,
            "trna_count_min" => trna_count_min,
            "n50_min" => n50_min,
            "strain_heterogeneity_max" => sh_max,
            "TreeCluster_method" => method,
            "TreeCluster_threshold" => threshold,
            "n_clusters" => n_clusters,
            "singletons" => singletons,
            "doubletons" => doubletons,
            "chao1" => chao1(clusters_cm, singletons, doubletons),
            "n_clustered_genomes" => n_clustered_genomes
        ),
        "accessions" => collect(zip(joined_df.accession, joined_df.cluster))
    )

    output_json = joinpath(dir, "selected_accessions.json")
    open(output_json, "w") do io
        JSON3.pretty(io, filtered_hq_accessions, JSON3.AlignmentContext(alignment=:Left, indent=2, level=0, offset=0))
    end
    @info "Results saved to: $output_json"
    exit()

    selected = combine(groupby(joined_df, :Cluster), :accession => first => :accession)
    count_selected = nrow(selected)
    @info "Selected representatives: $count_selected"

    selected_accessions = map(selected.accession) do acc
        m = match(r"GC[AF]_\d+\.\d+", acc)
        return m === nothing ? acc : m.match
    end

    if isnothing(output)
        output = joinpath(dir, "selected_ncbi_accessions.txt")
    end

    if isfile(output)
        error("Output file $output already exists. Use -o to specify another name.")
    end

    mkpath(dirname(output))
    open(output, "w") do io
        foreach(println, io, selected_accessions)
    end

    @info "Results saved to: $output"
    @info "Total accessions: $(length(selected_accessions))"
end

main()