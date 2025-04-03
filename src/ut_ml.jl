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
    class_counts = ClassWeights[]
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
        push!(class_counts, class_weights(chromosome))
    end
    return preallocated_chromosomes, class_counts
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
    
    labels_padded, class_counts = CDS_borders_both_strands(gff_data, chrom_names, chrom_lengths, side)
    labels_onehot_labels = (0,1,2,3)
    labels_encoded = onehotbatch.(labels_padded, Ref(labels_onehot_labels))
    
    return (dna_encoded, labels_encoded, class_counts)
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


const ChromosomeOneHot = OneHotArrays.OneHotMatrix{UInt32, Vector{UInt32}}
const LabelsOneHot = OneHotArrays.OneHotMatrix{UInt32, Vector{UInt32}}

const ClassWeights = Vector{Float32}
const GenomePrepared = Tuple{
    Vector{ChromosomeOneHot}, 
    Vector{LabelsOneHot},
    Vector{ClassWeights}
}
const IndexedGenome = Tuple{Int, GenomePrepared}


struct GenomeDataset
    genome_dirs::Vector{String}
    genome_counts::Int
    chromosome_counts::Vector{Int}
    pad::Int
    cds_side::Symbol
    max_cache_entries::Int

    _cache::Vector{IndexedGenome}
    _cum_counts::Vector{Int}

    function GenomeDataset(genome_dirs::Vector{T}; cds_side::Symbol=:starts, pad::Int=0, max_cache_entries::Int=10) where T <: AbstractString
        genome_counts = length(genome_dirs)
        chromosome_counts = count_chromosomes(genome_dirs)
        cum_counts = cumsum(chromosome_counts)

        new(
            string.(genome_dirs),
            genome_counts,
            chromosome_counts,
            pad,
            cds_side,
            max_cache_entries,
            IndexedGenome[],
            cum_counts
        )
    end
end

function Base.show(io::IO, ds::GenomeDataset)
    println(io, "GenomeDataset:")
    println(io, "  Genome Directories: ", first(ds.genome_dirs), ", ...")
    println(io, "  Genome Count      : ", ds.genome_counts)
    println(io, "  Total Chromosomes : ", length(ds))
    println(io, "  Padding           : ", ds.pad)
    println(io, "  CDS sides         : ", ds.cds_side)
    println(io, "  Cache Size        : ", length(ds._cache),"/",ds.max_cache_entries)
end

function _find_in_cache(cache::Vector{IndexedGenome}, idx::Int)
    for (genome_idx, genome) in cache
        if genome_idx == idx
            return genome
        end
    end
    return nothing
end

function _insert_cache!(ds::GenomeDataset, genome_index::Int, genome::GenomePrepared)
    length(ds._cache) >= ds.max_cache_entries && popfirst!(ds._cache)
    push!(ds._cache, (genome_index, genome))
end

function _get_or_load_cache(ds::GenomeDataset, idx::Int)
    cached = _find_in_cache(ds._cache, idx)
    !isnothing(cached) && return cached

    genome = process_genome_one_side(ds.genome_dirs[idx]; side=ds.cds_side, pad=ds.pad)
    _insert_cache!(ds, idx, genome)
    return genome
end

Base.length(ds::GenomeDataset) = last(ds._cum_counts)

function Base.getindex(ds::GenomeDataset, idx::Int)
    if idx < 1 || idx > length(ds)
        throw(BoundsError(ds, idx))
    end

    genome_idx = searchsortedfirst(ds._cum_counts, idx)
    prev_total = genome_idx == 1 ? 0 : ds._cum_counts[genome_idx - 1]
    chrom_idx = idx - prev_total

    genome = _get_or_load_cache(ds, genome_idx)
    return (genome[1][chrom_idx], genome[2][chrom_idx], genome[3][chrom_idx])
end

function Base.iterate(ds::GenomeDataset, state::Int)
    if state > length(ds)
        return nothing
    else
        return (ds[state], state + 1)
    end
end

function Base.iterate(ds::GenomeDataset)
    return Base.iterate(ds, 1)
end