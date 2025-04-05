using FastaIO
using Flux: onehotbatch, OneHotArrays

function CDS_borders_both_strands(gff_data, chrom_names, chrom_lengths, side::Symbol)
    @assert !isempty(gff_data) "Empty gff data"
    @assert !isempty(chrom_names) "Empty chromosomes list"
    @assert side in [:starts, :stops] "wrong side symbol: $side, must be either :starts or :stops"
    
    preallocated_chromosomes = [
        zeros(UInt8, chrom_length)
        for chrom_length in chrom_lengths
    ]
    for (chrom_name, chromosome) in zip(chrom_names, preallocated_chromosomes)
        cds_filter_pos = filter_gff_region(;
            sequence_header=chrom_name,
            regiontype="CDS",
            strand="+")
        cds_filter_neg = filter_gff_region(;
            sequence_header=chrom_name,
            regiontype="CDS",
            strand="-")

        CDS_regions_pos = cds_filter_pos(gff_data)
        CDS_regions_neg = cds_filter_neg(gff_data)

        if side == :starts
            ranges_from_GFF_records!(chromosome, CDS_regions_pos, :left; k=0x01)
            ranges_from_GFF_records!(chromosome, CDS_regions_neg, :right; k=0x02)
        else
            ranges_from_GFF_records!(chromosome, CDS_regions_pos, :right; k=0x01)
            ranges_from_GFF_records!(chromosome, CDS_regions_neg, :left; k=0x02)
        end
    end
    return preallocated_chromosomes
end

function class_weights(v::Vector{UInt8})
    L = length(v)
    n1 = count(x -> x == 0x01, v)
    n2 = count(x -> x == 0x02, v)
    n3 = count(x -> x == 0x03, v)
    return Float32[L-n1-n2-n3, n1, n2, n3] ./ L
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
    dna_onehot_labels = ('A', 'C', 'G', 'T', 'N')
    dna_encoded =  onehotbatch.(dna_padded, Ref(dna_onehot_labels), 'N')

    
    chrom_names = getindex.(fasta_data, 1) .|> split .|> first
    chrom_lengths = length.(fasta_sequences)
    gff_data = open_gff(gffname)
    
    labels_padded = CDS_borders_both_strands(gff_data, chrom_names, chrom_lengths, side)
    labels_onehot_labels = (0,1,2,3)
    labels_encoded = onehotbatch.(labels_padded, Ref(labels_onehot_labels))
    
    return (dna_encoded, labels_encoded)
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

function compute_confusion_matrices(true_labels::Vector{Int}, predictions::Vector{Int}, n_classes=4)
    if length(true_labels) != length(predictions)
        error("The length of true_labels and predictions must be the same.")
    end
    confusion_matrices = [zeros(Int, 2,2) for _ in 1:n_classes]
    for class in 1:n_classes
        TP = FP = FN = TN = 0
        for (t, p) in zip(true_labels, predictions)
            if t == class
                if p == class
                    TP += 1
                else
                    FN += 1
                end
            else
                if p == class
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