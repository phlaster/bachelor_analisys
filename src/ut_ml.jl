using FastaIO
using LuxCUDA
using Flux
using MLUtils
using Statistics
using ProgressMeter
using Functors
using Serialization
using StatsBase


function CDS_borders_both_strands(gff_data, chrom_names, chrom_lengths, side::Symbol)
    @assert !isempty(gff_data) "Empty gff data"
    @assert !isempty(chrom_names) "Empty chromosomes list"
    @assert side in [:starts, :stops] "wrong side symbol: $side, must be either :starts or :stops"
    
    preallocated_chromosomes = [
        zeros(Bool, chrom_length)
        for chrom_length in chrom_lengths
    ]
    for (chrom_name, chromosome) in zip(chrom_names, preallocated_chromosomes)
        cds_filter_pos = filter_gff_region(;
            sequence_header=chrom_name,
            regiontype="CDS",
            strand="+"
        )
        # cds_filter_neg = filter_gff_region(;
        #     sequence_header=chrom_name,
        #     regiontype="CDS",
        #     strand="-")

        CDS_regions_pos = cds_filter_pos(gff_data)
        # CDS_regions_neg = cds_filter_neg(gff_data)

        if side == :starts
            ranges_from_GFF_records!(chromosome, CDS_regions_pos, :left)
            # ranges_from_GFF_records!(chromosome, CDS_regions_neg, :right; k=0x02)
        else
            ranges_from_GFF_records!(chromosome, CDS_regions_pos, :right)
            # ranges_from_GFF_records!(chromosome, CDS_regions_neg, :left; k=0x02)
        end
    end
    return preallocated_chromosomes
end

function process_genome_one_side(genome_dir::T; side::Symbol=:starts, pad::Int=0) where T <: AbstractString
    @assert side in [:starts, :stops] "wrong side symbol: $side, must be either :starts or :stops"
    @assert isdir(genome_dir) "No such genomic directory"
    file_prefix = last(splitpath(genome_dir))
    fastaname = joinpath(genome_dir, "$file_prefix.fna")
    gffname = joinpath(genome_dir, "$file_prefix.gff")
    @assert all(isfile.([fastaname, gffname])) "Genomic files not found"

    
    fasta_data = readfasta(fastaname)
    fasta_sequences = getindex.(fasta_data, 2)
    dna_padded = add_pad.(collect.(uppercase.(fasta_sequences)), pad)
    dna_encoded =  Flux.onehotbatch.(dna_padded, Ref(('A', 'C', 'G', 'T', 'N')), 'N')

    
    chrom_names = getindex.(fasta_data, 1) .|> split .|> first
    chrom_lengths = length.(fasta_sequences)
    gff_data = open_gff(gffname)
    
    labels_padded = CDS_borders_both_strands(gff_data, chrom_names, chrom_lengths, side)
    # labels_encoded = Flux.onehotbatch.(labels_padded, Ref((0,1)))
    
    return (dna_encoded, labels_padded)
end

function count_chromosomes(genome_dir::T) where T <: AbstractString
    @assert isdir(genome_dir) "No such genomic directory"
    countfilename = only(filter(startswith("n_chroms="), readdir(genome_dir)))
    n_chroms = parse(Int, countfilename[10:end])
    return n_chroms
end

function count_chromosomes(genome_dirs::Vector{T}) where T <: AbstractString
    tasks = [Threads.@spawn count_chromosomes(dir) for dir in genome_dirs]
    results = fetch.(tasks)
    return results
end

function add_pad(v::Vector, pad::Int)
    pad < 1 && return v
    L = length(v)
    pad <= L && return vcat(v[end-pad+1:end], v, v[1:pad])

    leading = [v[mod1(i, L)] for i in L-pad+1:L]
    trailing = [v[mod1(i, L)] for i in 1:pad]
    
    return vcat(leading, v, trailing)
end

_remove_singular_dims(x) = dropdims(x, dims=(2,3))
function create_model(; window_size::Int, input_channels::Int=5)
    Chain(
        Conv((window_size,), input_channels => 128, pad=0, bias=false),
        BatchNorm(128, relu),
        Dropout(0.15),
        
        Conv((5,), 128 => 64, pad=2, bias=false),
        BatchNorm(64, relu),
        Dropout(0.15),

        Conv((3,), 64 => 64, pad=1, bias=false),
        BatchNorm(64, relu),
        
        Conv((1,), 64 => 1),
        _remove_singular_dims,
        sigmoid
    )
end


function split_into_chunks(seq, labels, max_chunk_size)
    seq_length = size(seq, 2)
    labels_length = length(labels)
    pad = (seq_length - labels_length) ÷ 2
    seq_length <= max_chunk_size && return [(seq, labels)]

    nchunks = labels_length ÷ max_chunk_size + 1
    chunk_lengths_approx = ceil(Int, labels_length/nchunks)
    chunk_coords = Base.Iterators.partition(1:labels_length, chunk_lengths_approx)

    res = @views (
        (seq[:, (coords.start):(coords.stop+2pad)], labels[coords])
        for coords in chunk_coords
    )
   return res 
end

_transform_X(seq) = reshape(permutedims(seq, (2, 1)), size(seq, 2), size(seq, 1), 1)
_transform_y(labels) = permutedims(labels, (2, 1))
_add_grads(a, b) = isnothing(a) ? b : isnothing(b) ? a : a .+ b
_scale_grads(x, factor) = isnothing(x) ? nothing : x .* factor

function _train_epoch!(model, dataset, opt, loss_function, max_chunk_size, dev, chunk_skip_coeff)
    @info "Optim: $(opt.layers[1].weight.rule)"
    
    batch_losses = Float32[]
    n_labels = 0
    n_positive = 0
    n_chunks = 0
    skipped_chunks = 0
    @showprogress desc="Epoch progress:" barlen=50 color=:blue dt=1 showspeed=true for (seq, labels) in dataset
        accumulated_grads = nothing
        accum_counter = 0

        for (chunk_X, chunk_y) in split_into_chunks(seq, labels, max_chunk_size)
            n_chunks += 1
            n_labels_chunk = length(labels)
            n_chunk_positive = sum(labels)
            n_labels += n_labels_chunk
            n_positive += n_chunk_positive
            if n_chunk_positive/n_labels_chunk < n_positive/n_labels * chunk_skip_coeff
                # comparing per-chunk positive class frequency (PCF) with overall PCF
                # skipping chunks with low PCF
                skipped_chunks += 1
                continue
            end

            X = _transform_X(chunk_X) |> dev
            y = chunk_y               |> dev

            loss, grads = Flux.withgradient(model) do m
                ŷ = m(X)
                loss_function(ŷ, y)
            end
            accumulated_grads = Functors.fmap(_add_grads, accumulated_grads, grads[1])
            accum_counter += 1
            push!(batch_losses, loss)
        end

        if accum_counter > 0
            scaled_grads = Functors.fmap(x -> _scale_grads(x, inv(accum_counter)), accumulated_grads)
            Flux.update!(opt, model, scaled_grads)
        end
    end
    epoch_loss = mean(batch_losses)
    @info "Chunks skipped: $skipped_chunks/$n_chunks = $(round(skipped_chunks/n_chunks, digits=3))"
    @info "Mean loss     : $epoch_loss"
    return epoch_loss, skipped_chunks, n_chunks
end

function train_model(model, ds_train::GenomeDataset, ds_test::GenomeDataset;
    epochs::Int=1,
    lr::Float64=0.005,
    dev=gpu,
    floss_gamma::Float64=3.5,
    loss_lambda::Float64=0.1,
    decay_factor::Float64=0.6,
    max_chunk_size::Int=2*10^5,
    chunk_skip_coeff::Float64=0.0,
    savedir=".")
    
    @assert all([
        epochs>0,
        0<lr≤1,
        0≤floss_gamma,
        0≤loss_lambda,
        0<decay_factor≤1,
        0<max_chunk_size≤10^6,
        0≤chunk_skip_coeff
    ]) "Wrong parameter value"
    @assert isdir(savedir) "Wrong directory name for model saving"

    model = model |> dev

    dev_label = CUDA.device(first(model.layers).weight) |> string

    opt = Flux.setup(Adam(lr), model)

    function loss_function(ŷ, y; gamma=floss_gamma, min_distance=60, lambda=loss_lambda)
        focal_loss = Flux.Losses.binary_focal_loss(ŷ, y; gamma=gamma)
        predictions = findall(ŷ .> 0.5)
        distances = diff(predictions)
        spatial_penalty = sum(exp.(-distances / min_distance))
    
        total_loss = focal_loss + lambda * spatial_penalty
        return total_loss
    end

    losses = Float32[]
    lrs = Float64[]
    metrics = Dict{Int64, NamedTuple}[]
    cms = Matrix{Int}[]
    
    current_lr = lr

    for epoch in 1:epochs
        @info "Epoch $epoch/$epochs"
        dumpname = joinpath(savedir, "$(dev_label)_epoch_$(lpad(epoch, 3, '0')).flux")
        
        epoch_loss, skipped_chunks, n_chunks = _train_epoch!(
            model, ds_train, opt, loss_function, max_chunk_size, dev, chunk_skip_coeff
        )
        (class_metrics, conf_mtr, classes), fp_shifts, gt_s2s_ranges, pred_s2s_ranges = evaluate_bin_class_model(
            model, ds_test; dev=dev, max_chunk_size=max_chunk_size
        )
        
        push!(metrics, class_metrics)
        push!(cms, conf_mtr)
        push!(losses, epoch_loss)
        push!(lrs, current_lr)
        
        dump_data = (
            # Frozen state
            model=model |> cpu,
            opt=opt |> cpu,

            # Single-value vars
            epoch=epoch,
            classes=classes,
            pad=ds_train.pad,
            device=dev_label,
            gamma=floss_gamma,
            lambda=loss_lambda,
            decay=decay_factor,
            chunk=max_chunk_size,
            chunk_skip=chunk_skip_coeff,
            n_chunks=n_chunks,
            skipped_chunks=skipped_chunks,
            dir=savedir,
            
            # Train history vars
            lr=lrs,
            metrics=metrics,
            cm=cms,
            loss=losses,

            # Aggregated countmaps
            fp_shifts=fp_shifts,
            cds_ranges_true=gt_s2s_ranges,
            cds_ranges_model=pred_s2s_ranges
        )
        @info "Prec  : $(round(class_metrics[1].prec, digits=4))"
        @info "Recall: $(round(class_metrics[1].rec, digits=4))"
        @info "F1    : $(round(class_metrics[1].f1, digits=4))"
        mkpath(dirname(dumpname))
        serialize(dumpname, dump_data)
        @info "Saved training dump: $dumpname"

        current_lr = lr * decay_factor^epoch
        Flux.adjust!(opt, current_lr)
    end
end

function dotrain_model(model, ds_train::GenomeDataset, ds_test::GenomeDataset;
    epochs::Int,
    dev,
    opt,
    elapsed_epochs::Int,
    floss_gamma::Float64,
    loss_lambda::Float64,
    decay_factor::Float64,
    max_chunk_size::Int,
    chunk_skip_coeff::Float64,
    losses::Vector{Float32},
    lrs::Vector{Float64},
    metrics::Vector{Dict{Int64, NamedTuple}},
    cms::Vector{Matrix{Int}},
    savedir)   
    @assert all([
        epochs>0,
        length(losses)==length(lrs)==length(metrics)==length(cms)==elapsed_epochs,
        0 ≤ floss_gamma,
        0 ≤ loss_lambda,
        0 < decay_factor ≤1,
        0 < max_chunk_size ≤10^6,
        0 ≤ chunk_skip_coeff
    ]) "Wrong parameter value"

    model = model |> dev
    opt = opt |> dev

    dev_label = CUDA.device(first(model.layers).weight) |> string
    function loss_function(ŷ, y; gamma=floss_gamma, min_distance=60, lambda=loss_lambda)
        focal_loss = Flux.Losses.binary_focal_loss(ŷ, y; gamma=gamma)
        predictions = findall(ŷ .> 0.5)
        distances = diff(predictions)
        spatial_penalty = sum(exp.(-distances / min_distance))
    
        total_loss = focal_loss + lambda * spatial_penalty
        return total_loss
    end
    current_lr = last(lrs) * decay_factor

    start_epoch = elapsed_epochs+1
    last_epoch = start_epoch + epochs-1
    for epoch in start_epoch:last_epoch
        @info "Epoch $epoch/$last_epoch"
        dumpname = joinpath(savedir, "$(dev_label)_epoch_$(lpad(epoch, 3, '0')).flux")
        
        epoch_loss, skipped_chunks, n_chunks = _train_epoch!(
            model, ds_train, opt, loss_function, max_chunk_size, dev, chunk_skip_coeff
        )
        (class_metrics, conf_mtr, classes), fp_shifts, gt_s2s_ranges, pred_s2s_ranges = evaluate_bin_class_model(
            model, ds_test; dev=dev, max_chunk_size=max_chunk_size
        )
        
        push!(metrics, class_metrics)
        push!(cms, conf_mtr)
        push!(losses, epoch_loss)
        push!(lrs, current_lr)
        
        dump_data = (
            # Frozen state
            model=model |> cpu,
            opt=opt |> cpu,

            # Single-value vars
            epoch=epoch,
            classes=classes,
            pad=ds_train.pad,
            device=dev_label,
            gamma=floss_gamma,
            lambda=loss_lambda,
            decay=decay_factor,
            chunk=max_chunk_size,
            chunk_skip=chunk_skip_coeff,
            n_chunks=n_chunks,
            skipped_chunks=skipped_chunks,
            dir=savedir,
            
            # Train history vars
            lr=lrs,
            metrics=metrics,
            cm=cms,
            loss=losses,

            # Aggregated countmaps
            fp_shifts=fp_shifts,
            cds_ranges_true=gt_s2s_ranges,
            cds_ranges_model=pred_s2s_ranges
        )
        @info "Prec  : $(round(class_metrics[1].prec, digits=4))"
        @info "Recall: $(round(class_metrics[1].rec, digits=4))"
        @info "F1    : $(round(class_metrics[1].f1, digits=4))"
        mkpath(dirname(dumpname))
        serialize(dumpname, dump_data)
        @info "Saved training dump: $dumpname"

        current_lr *= decay_factor
        Flux.adjust!(opt, current_lr)
    end
end

function true_distances_inside_chunk(v)
    if eltype(v) != Bool
        throw(TypeError(:true_distances_inside_chunk, Bool, eltype(v)))
    end
    v |> findall |> diff
end

function false_positive_stats(ground_truth, predicted)
    @assert length(ground_truth) == length(predicted) "Vectors must be of equal length"
    
    ground_truth_idx = findall(ground_truth)    
    matched_ground_truth = Set{Int}()
    fp_candidates = Int[]
    
    for i in eachindex(ground_truth)
        if predicted[i]
            if ground_truth[i]
                push!(matched_ground_truth, i)
            else
                push!(fp_candidates, i)
            end
        end
    end

    unmatched_ground_truth = [i for i in ground_truth_idx if i ∉ matched_ground_truth]
    matched_distances = Int[]
   
    candidates = NTuple{4, Int}[]
    for fp in fp_candidates
        for gt in unmatched_ground_truth
            signed_dist = gt - fp   # signed distance (negative if gt < fp)
            abs_dist = abs(signed_dist)
            push!(candidates, (abs_dist, fp, gt, signed_dist))
        end
    end
    sorted_candidates = sort(candidates, by = first)
    
    unmatched_set = Set(unmatched_ground_truth)
    matched_fp = Set{Int}()

    for candidate in sorted_candidates
        _, fp, gt, signed_dist = candidate
        if (fp ∉ matched_fp) && (gt ∈ unmatched_set)
            # Check if there's any matched TP (from matched_ground_truth) between fp and gt
            has_intermediate = false
            min_pos = min(fp, gt)
            max_pos = max(fp, gt)
            for m_gt in matched_ground_truth
                if m_gt > min_pos && m_gt < max_pos
                    has_intermediate = true
                    break
                end
            end
            if !has_intermediate
                push!(matched_distances, signed_dist)
                push!(matched_fp, fp)
                delete!(unmatched_set, gt)
            end
        end
    end
    return matched_distances
end

function evaluate_bin_class_model(model, dataset; dev=gpu, max_chunk_size=2*10^5)
    function addvals!(d1::Dict{T, Int64}, d2::Dict{T, Int64}) where T
        for k2 in keys(d2)
            if haskey(d1, k2)
                d1[k2] += d2[k2]
            else
                d1[k2] = d2[k2]
            end
        end
    end
    
    cm_classes = zeros(Int, 2,2)
    fp_shifts = Dict{Int, Int}()
    gt_s2s_ranges = Dict{Int, Int}() # start-to-start
    pred_s2s_ranges = Dict{Int, Int}()
    
    classes = (0,1)
    @showprogress desc="Evaluating:" barlen=50 color=:white showspeed=true for (seq, labels) in dataset
        for (seq_chunk, labels_chunk) in split_into_chunks(seq, labels, max_chunk_size)        
            X = _transform_X(seq_chunk) |> dev
            y_pred = model(X) .|> >=(0.5f0) |> cpu
            cm_classes .+= multiclass_confusion(labels_chunk, y_pred; classes=classes)[1]

            addvals!(fp_shifts,       false_positive_stats(labels_chunk, y_pred)|> countmap)
            addvals!(gt_s2s_ranges,   true_distances_inside_chunk(labels_chunk) |> countmap)
            addvals!(pred_s2s_ranges, true_distances_inside_chunk(y_pred)       |> countmap)
        end
    end
    return metrics_from_cm(cm_classes, classes), fp_shifts, gt_s2s_ranges, pred_s2s_ranges
end

