

function create_sankey(steps::Vector{Pair{Int64, String}}; savefile=nothong)
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
    if !isnothing(savefile)
        PlotlyJS.savefig(plt, savefile)
        @info "Sankey diagrame saved to $savefile"
    end
    return plt
end
