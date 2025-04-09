using FastaIO
using LuxCUDA
using Flux
using MLUtils
using Statistics
using ProgressMeter
using Functors


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
    tasks = [@async count_chromosomes(dir) for dir in genome_dirs]
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

function compute_confusion_matrices(true_labels::Vector{T}, predictions::Vector{T}, n_classes=2) where T <: Integer
    if length(true_labels) != length(predictions)
        error("The length of true_labels and predictions must be the same.")
    end
    confusion_matrices = [zeros(Int, 2,2) for _ in 1:n_classes]
    for class in 1:n_classes
        TP = FP = FN = TN = 0
        for (t, p) in zip(true_labels, predictions)
            if t+1 == class
                if p+1 == class
                    TP += 1
                else
                    FN += 1
                end
            else
                if p+1 == class
                    FP += 1
                else
                    TN += 1
                end
            end
        end        
        confusion_matrices[class] += [
            TP FP
            FN TN
        ]
    end
    return confusion_matrices
end

function compute_metrics(conf_mat::Matrix{Int})
    TP = conf_mat[1, 1]
    FP = conf_mat[1, 2]
    TN = conf_mat[2, 2]
    FN = conf_mat[2, 1]

    cm = ConfusionMTR("Confusion", (TP=TP,FP=FP,TN=TN,FN=FN))
    return cm
end


function create_model(; window_size::Int, input_channels::Int=5)
    _remove_singular_dims(x) = dropdims(x, dims=(2,3))
    
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

function train_model!(model, dataset::GenomeDataset;
    epochs=1:1,
    lr=0.001,
    device=gpu,
    max_chunk_size=2*10^5
)
    function _add_grads(a, b)
        if isnothing(a)
            if isnothing(b)
                return nothing
            else
                return b
            end
        elseif isnothing(b)
            return a
        end
        return a .+ b
    end
    _scale_grads(x, factor) = isnothing(x) ? nothing : x .* factor
    
    opt_state = Flux.setup(Adam(lr), model)
    epoch_losses = Float32[]
    epoch_loss = Inf

    for epoch in epochs
        batch_losses = Float32[]
        @showprogress desc="Epoch $epoch/$(epochs.stop)" for (seq, labels) in dataset
            accumulated_grads = nothing
            accum_counter = 0

            for (seq_chunk, labels_chunk) in split_into_chunks(seq, labels, max_chunk_size)
                X = _transform_X(seq_chunk) |> device
                y = labels_chunk            |> device

                loss, grads = Flux.withgradient(model) do m
                    ŷ = m(X)
                    Flux.Losses.binary_focal_loss(ŷ, y; gamma=3)
                end

                accumulated_grads = Functors.fmap(_add_grads, accumulated_grads, grads[1])
                accum_counter += 1
                push!(batch_losses, loss)
            end

            if accum_counter > 0
                scaled_grads = Functors.fmap(x -> _scale_grads(x, inv(accum_counter)), accumulated_grads)
                Flux.update!(opt_state, model, scaled_grads)
            end
        end

        epoch_loss = mean(batch_losses)
        
        @info "Mean loss for epoch $epoch: $epoch_loss"
        push!(epoch_losses, epoch_loss)
    end
    return epoch_losses
end

function evaluate_model(model, dataset; device=gpu, max_chunk_size=2*10^5, n_classes=2)

    cm_classes = [zeros(Int, 2,2) for _ in 1:n_classes]

    @showprogress desc = "Evaluating" for (seq, labels) in dataset
        for (seq_chunk, labels_chunk) in split_into_chunks(seq, labels, max_chunk_size)        
            X = _transform_X(seq_chunk) |> device
            y_pred = model(X) .|> >=(0.5f0) |> cpu
            cm_classes .+= compute_confusion_matrices(collect(labels_chunk), y_pred)
        end
    end

    return compute_metrics.(cm_classes)
end