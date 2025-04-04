const ChromosomeOneHot = OneHotArrays.OneHotMatrix{UInt32, Vector{UInt32}}
const LabelsOneHot = OneHotArrays.OneHotMatrix{UInt32, Vector{UInt32}}

const GenomePrepared = Tuple{
    Vector{ChromosomeOneHot}, 
    Vector{LabelsOneHot},
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
    return (genome[1][chrom_idx], genome[2][chrom_idx])
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