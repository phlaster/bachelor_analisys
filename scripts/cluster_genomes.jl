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
using Chain
using PlotlyJS

const METADATA_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz"
const TREE_URL = "https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_r220.tree.gz"
const GENBANK_PATTERN = r"GC[AF]_\d+\.\d+"
extract_accession(acc) = match(GENBANK_PATTERN, acc).match


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--InputDirectory", "-I"
            help = "Output directory for tree and metadata files (required)"
            required = true
        "--OutputDirectory", "-O"
            help = "Output file name for selected accessions"
        "--Accessions", "-A"
            help = "Flag to include accessions into resulting JSON"
            action = :store_true

        "--TreeClusterTreshold", "-T"
            help = "Clustering threshold for TreeCluster.py (required)"
            default = 0.2463
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
            default = 98.0
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
            default = 40
            arg_type = Int64
            range_tester = x->x≥0
        "--N50Min", "-N"
            help = "N50 score for a genome"
            default = 50_000
            arg_type = Int64
            range_tester = x->x≥1
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

function create_sankey(steps::Vector{Pair{Int64, String}})
    # Number of steps in the data pipeline.
    n = length(steps)
    # We want a kept node for every step plus one extra "Out" node.
    n_keep = n + 1

    # These arrays will hold all node properties.
    node_labels = String[]
    node_colors = String[]
    node_x = Float64[]
    node_y = Float64[]

    # Arrays for link properties.
    link_source = Int[]
    link_target = Int[]
    link_value = Int[]

    # We will interleave kept nodes and, if needed, a discard node right before a kept node.
    # To position the kept nodes in a horizontal row we set:
    kept_x(k) = (k - 1) / (n_keep - 1)  # k = 1,2,…,n_keep
    kept_y = 0.45

    # For discard nodes we want them “branched” off the kept node that follows.
    # Here we choose a small x offset and a different y so that they don’t overlap.
    discard_x(k) = kept_x(k) - 0.06  # associated with kept node order k (k ≥ 2)
    discard_y = 0.8

    # We'll also record the final node index (PlotlyJS sankey uses 0-based indexing)
    # and store the indices for kept nodes.
    final_node_index = 0
    kept_indices = Int[]

    # Add the very first kept node ("Before filtering"); its color is blue.
    push!(node_labels, steps[1].second)
    push!(node_colors, "#0000ff90")
    push!(node_x, kept_x(1))
    push!(node_y, kept_y)
    push!(kept_indices, final_node_index)
    final_node_index += 1

    # For each filter stage (from step 2 up to the last step)
    for i in 2:n
        # Calculate how many samples were discarded at this stage.
        discarded = steps[i-1].first - steps[i].first
        if discarded > 0
            # Compute percentage (with two decimals)
            pct = 100 * discarded / steps[i-1].first
            label = "-$discarded (-$(round(pct, digits=2))%)"
            # Add a discard node (color red) associated with the next kept node (order = i)
            push!(node_labels, label)
            push!(node_colors, "#ff000070")
            push!(node_x, discard_x(i))
            push!(node_y, discard_y)
            discard_node_index = final_node_index
            final_node_index += 1

            # Create a link from the previous kept node (i–1) to this discard node.
            push!(link_source, kept_indices[end])
            push!(link_target, discard_node_index)
            push!(link_value, discarded)
        end
        # Now add the kept node for this filter stage (color blue).
        push!(node_labels, steps[i].second)
        push!(node_colors, "#0000ff90")
        push!(node_x, kept_x(i))
        push!(node_y, kept_y)
        push!(kept_indices, final_node_index)
        # Create the main link from the previous kept node to this one.
        push!(link_source, kept_indices[end-1])
        push!(link_target, final_node_index)
        push!(link_value, steps[i].first)
        final_node_index += 1
    end
    push!(node_labels, "Filtered")
    push!(node_colors, "#00cc0099")
    push!(node_x, kept_x(n_keep))
    push!(node_y, kept_y)
    push!(kept_indices, final_node_index)
    # Create the last link from the last filter stage to "Out".
    push!(link_source, kept_indices[end-1])
    push!(link_target, final_node_index)
    push!(link_value, steps[end].first)
    final_node_index += 1

    # Build the sankey diagram.
    trace = sankey(
        node = attr(
            label = node_labels,
            color = node_colors,
            x = node_x,
            y = node_y,
            pad = 40,
            thickness = 10,
            line = attr(color = "black", width = 0.7),
        ),
        link = attr(
            source = link_source,
            target = link_target,
            value = link_value
        )
    )
	
    plt = plot(trace,
		Layout(title_text = "Applied filters", font_size = 15, margin = attr(b = 100, r=100, l=100, t=100))
	)
    return plt
end

function pick_one_per_cluster(df::DataFrame; skip=String[])
    df_filtered = filter(row -> !(row.accession in skip), df)
    groups = groupby(df_filtered, :cluster)
    picks = [first(g.accession) for g in groups]
    return picks
end


function main()
    check_dependencies()
    args = parse_commandline()

    dir_in = abspath(args["InputDirectory"])
    dir_out = isnothing(args["OutputDirectory"]) ? dir_in : args["OutputDirectory"]
    include_accessions = args["Accessions"]

    threshold = args["TreeClusterTreshold"]
    method = args["TreeClusterMethod"]

    completeness_min = args["CompletenessMin"]
    contamination_max = args["ContaminationMax"]
    contig_count_max = args["ContigCountMax"]
    trna_count_min = args["tRNACountMin"]
    n50_min = args["N50Min"]
    
    mkpath(dir_in)
    mkpath(dir_out)
    
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
    output_sankey_filter = joinpath(dir_out, "sankey_filter.html")
    open(output_sankey_filter, "w") do io
	    PlotlyBase.to_html(io, create_sankey(filter_steps).plot)
	end
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
        "parameters" => Dict(
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
        "accessions" => include_accessions ? collect(zip(joined_df.accession, joined_df.cluster)) : []
    )

    output_json = joinpath(dir_out, "filtering_report.json")
    open(output_json, "w") do io
        JSON3.pretty(io, filtered_hq_accessions, JSON3.AlignmentContext(alignment=:Left, indent=2, level=0, offset=0))
    end
    @info "Report saved to: $output_json"

    output_all_accessions = joinpath(dir_out, "filtered_accessions.txt")
    open(output_all_accessions, "w") do io
        for acc in joined_df.accession
            println(io, acc)
        end
    end
    @info "All filtered accession saved to $output_all_accessions"

    
    train_accs = pick_one_per_cluster(joined_df)
    validation_accs = pick_one_per_cluster(joined_df; skip=train_accs)

    train_filename = joinpath(dir_out, "train.txt")
    validation_filename = joinpath(dir_out, "validation.txt")

    open(train_filename, "w") do io
        foreach(acc -> println(io, acc), train_accs)
    end
    @info "Train accessions written to $train_filename"

    open(validation_filename, "w") do io
        foreach(acc -> println(io, acc), validation_accs)
    end
    @info "Validation accessions written to $validation_filename"
end

main()