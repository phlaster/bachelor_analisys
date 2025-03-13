using PlotlyJS
using ColorSchemes

function create_sankey(steps::Vector{Pair{Int64, String}}; savefile=nothong)
    n = length(steps)
    n_keep = n + 1

    node_labels = String[]
    node_colors = String[]
    node_x = Float64[]
    node_y = Float64[]

    link_source = Int[]
    link_target = Int[]
    link_value = Int[]

    kept_x(k) = (k - 1) / (n_keep - 1)
    kept_y = 0.45

    discard_x(k) = kept_x(k) - 0.06
    discard_y = 0.8

    final_node_index = 0
    kept_indices = Int[]

    push!(node_labels, steps[1].second)
    push!(node_colors, "#0000ff90")
    push!(node_x, kept_x(1))
    push!(node_y, kept_y)
    push!(kept_indices, final_node_index)
    final_node_index += 1

    for i in 2:n
        discarded = steps[i-1].first - steps[i].first
        if discarded > 0
            pct = 100 * discarded / steps[i-1].first
            label = "-$discarded (-$(round(pct, digits=2))%)"
            push!(node_labels, label)
            push!(node_colors, "#ff000070")
            push!(node_x, discard_x(i))
            push!(node_y, discard_y)
            discard_node_index = final_node_index
            final_node_index += 1

            push!(link_source, kept_indices[end])
            push!(link_target, discard_node_index)
            push!(link_value, discarded)
        end
        push!(node_labels, steps[i].second)
        push!(node_colors, "#0000ff90")
        push!(node_x, kept_x(i))
        push!(node_y, kept_y)
        push!(kept_indices, final_node_index)
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
    push!(link_source, kept_indices[end-1])
    push!(link_target, final_node_index)
    push!(link_value, steps[end].first)
    final_node_index += 1

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
    if !isnothing(savefile)
        PlotlyJS.savefig(plt, savefile)
        @info "Sankey diagram saved to $savefile"
    end
    return plt
end

function summaryplots(main_df; train_df, savename)
    # Calculate axis limits with 10% padding for scatter plots
    genome_min, genome_max = extrema(main_df.genome_size)
    gene_min, gene_max = extrema(main_df.gene_count)
    genome_padding = 0.1 * (genome_max - genome_min)
    gene_padding = 0.1 * (gene_max - gene_min)
    
    xlims = (genome_min - genome_padding, genome_max + genome_padding)
    ylims = (gene_min - gene_padding, gene_max + gene_padding)

    # Existing phylum processing and color mapping code remains the same
    phylum_counts = combine(groupby(main_df, :phylum), nrow => :count)
    sort!(phylum_counts, :count, rev=true)
    total = sum(phylum_counts.count)
    cumsum_counts = cumsum(phylum_counts.count)
    threshold = 0.9 * total
    cutoff_idx = findfirst(≥(threshold), cumsum_counts)
    
    major_phylums = if isnothing(cutoff_idx)
        unique(main_df.phylum)
    else
        phylum_counts.phylum[1:cutoff_idx]
    end

    color_scheme = ColorSchemes.flag_usmi
    color_map = Dict(
        phylum => get(color_scheme, i / length(major_phylums))
        for (i, phylum) in enumerate(major_phylums)
    )

    # Existing create_traces function remains the same
    function create_traces(df, name_prefix)
        traces = []
        for phylum in major_phylums
            sub_df = filter(row -> row.phylum == phylum, df)
            p_count = nrow(sub_df)
            p_percent = round(100nrow(sub_df)/nrow(df), digits=2)
            color = color_map[phylum]
            tr = scatter(x=sub_df.genome_size, y=sub_df.gene_count,
                mode="markers", name="$p_count($p_percent%) $phylum",
                marker=attr(
                    color="rgb($(color.r), $(color.g), $(color.b))",
                    opacity = 0.5,
                ),
                legendgroup=phylum,
                showlegend=true,
            )
            push!(traces, tr)
        end
        
        other_df = filter(row -> !(row.phylum in major_phylums), df)
        if nrow(other_df) > 0
            otr_count = nrow(other_df)
            otr_percent = round(100nrow(other_df)/nrow(df), digits=2)
            tr = scatter(x=other_df.genome_size, y=other_df.gene_count,
                mode="markers", name="$otr_count($otr_percent%) Other",
                marker=attr(
                    color="#CCCCCC",
                    opacity = 0.45,
                ),
                legendgroup="Other",
                showlegend=true
            )
            push!(traces, tr)
        end
        return traces
    end

    # Existing trace creation remains the same
    main_traces = create_traces(main_df, "Main")
    train_traces = isnothing(train_df) ? 
        [scatter(x=[NaN], y=[NaN], mode="markers", name="No Training Data")] :
        create_traces(train_df, "Train")

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[
            "Gene Count Histogram" "Genome Length vs Gene Count (Main), $(nrow(main_df))"
            "GC Percentage Histogram" "Genome Length vs Gene Count (Train), $(nrow(train_df))"
        ]
    )

    # Add traces with axis labels
    add_trace!(fig, histogram(x=main_df.gene_count, showlegend=false), row=1, col=1)
    add_trace!(fig, histogram(x=main_df.gc_percentage, showlegend=false), row=1, col=2)

    for trace in main_traces
        add_trace!(fig, trace, row=2, col=1)
    end
    for trace in train_traces
        add_trace!(fig, trace, row=2, col=2)
    end

    # Set axis labels and limits
    relayout!(fig, 
        title_text="Dataset and train subset statistics",
        legend=attr(x=1.05, y=0.5),
        margin=attr(r=200),
        xaxis1_title="Gene Count",
        yaxis1_title="Frequency",
        xaxis2_title="GC Percentage",
        yaxis2_title="Frequency",
        xaxis3_title="Genome Length",
        yaxis3_title="Gene Count",
        xaxis4_title="Genome Length",
        yaxis4_title="Gene Count",
        xaxis3_range=xlims,
        yaxis3_range=ylims,
        xaxis4_range=xlims,
        yaxis4_range=ylims
    )

    if !isnothing(savename)
        PlotlyJS.savefig(fig, savename)
        @info "Plots saved to $savename"
    end
    return fig
end