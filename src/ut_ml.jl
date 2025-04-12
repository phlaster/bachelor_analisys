using FastaIO
using LuxCUDA
using Flux
using MLUtils
using Statistics
using ProgressMeter
using Functors
using Serialization


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
        Dropout(0.1),
        
        Conv((3,), 128 => 64, pad=1, bias=false),
        BatchNorm(64, relu),
        Dropout(0.1),

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
                # skipping chunks with low positive class representation
                # comparing chunk positive class frequency (pcf) with overall pcf
                # logs are used to avoid integer overflow
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
    floss_gamma::Float64=3.0,
    decay_factor::Float64=0.6,
    max_chunk_size::Int=2*10^5,
    chunk_skip_coeff::Float64=0.0,
    savedir=".",
)   
    @assert all([
        epochs>0,
        0<lr≤1,
        0≤floss_gamma,
        0<decay_factor≤1,
        0<max_chunk_size≤10^6,
        0≤chunk_skip_coeff
    ]) "Wrong parameter value"
    @assert isdir(savedir) "Wrong directory name for model saving"

    model = model |> dev

    dev_label = CUDA.device(first(model.layers).weight) |> string

    opt = Flux.setup(Adam(lr), model)
    loss_function(ŷ, y) = Flux.Losses.binary_focal_loss(ŷ, y; gamma=floss_gamma)

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
        class_metrics, conf_mtr, classes = evaluate_bin_class_model(model, ds_test; dev=dev, max_chunk_size=max_chunk_size)
        
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
            loss=losses
        )
        @info "Prec  : $(round(class_metrics[1].prec, digits=4))"
        @info "Recall: $(round(class_metrics[1].rec, digits=4))"
        @info "F1    : $(round(class_metrics[1].f1, digits=4))"
        
        serialize(dumpname, dump_data)
        @info "Saved training dump: $dumpname"

        current_lr = lr * decay_factor^epoch
        Flux.adjust!(opt, current_lr)
    end
end

function evaluate_bin_class_model(model, dataset; dev=gpu, max_chunk_size=2*10^5)
    cm_classes = zeros(Int, 2,2)
    classes = (0,1)
    @showprogress desc="Evaluating:" barlen=50 color=:white showspeed=true for (seq, labels) in dataset
        for (seq_chunk, labels_chunk) in split_into_chunks(seq, labels, max_chunk_size)        
            X = _transform_X(seq_chunk) |> dev
            y_pred = model(X) .|> >=(0.5f0) |> cpu
            cm_classes .+= multiclass_confusion(labels_chunk, y_pred; classes=classes)[1]
        end
    end
    return metrics_from_cm(cm_classes, classes)
end