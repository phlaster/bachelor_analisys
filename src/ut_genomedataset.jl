using Flux: OneHotArrays

const Chromosome_X = OneHotArrays.OneHotMatrix{UInt32, Vector{UInt32}}
const Labels_y = Vector{Bool}

const GenomeSampleTuple = Tuple{
    Vector{Chromosome_X}, 
    Vector{Labels_y},
}

struct GenomeDataset
    genome_dirs::Vector{String}
    genome_counts::Int
    chromosome_counts::Vector{Int}
    pad::Int
    cds_side::Symbol
    
    _cache::Vector{Union{Nothing, GenomeSampleTuple}}
    _cum_counts::Vector{Int}

    function GenomeDataset(genome_dirs::Vector{String};
        cds_side::Symbol=:starts,
        pad::Int=0
    )
    
        genome_counts = length(genome_dirs)
        chromosome_counts = count_chromosomes(genome_dirs)
        cum_counts = cumsum(chromosome_counts)
        cache = Vector{Union{Nothing, GenomeSampleTuple}}(nothing, genome_counts)
        
        Threads.@spawn begin
            for idx in 1:genome_counts
                cache[idx] = process_genome_one_side(
                    genome_dirs[idx];
                    side=cds_side,
                    pad=pad
                )
            end
        end
        
        new(
            genome_dirs,
            genome_counts,
            chromosome_counts,
            pad,
            cds_side,
            cache,
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
    last_cached = findfirst(isnothing, ds._cache)
    last_cached = isnothing(last_cached) ? length(ds._cache) : last_cached-1
    println(io, "  Cached            : ", round(Int, 100*last_cached/length(ds._cache)), "%")
end

# function _find_in_cache(cache::Vector{IndexedGenome}, idx::Int)
#     for (genome_idx, genome) in cache
#         if genome_idx == idx
#             return genome
#         end
#     end
#     return nothing
# end

# function _insert_cache!(ds::GenomeDataset, genome_index::Int, genome::GenomeSampleTuple)
#     length(ds._cache) >= ds.max_cache_entries && popfirst!(ds._cache)
#     push!(ds._cache, (genome_index, genome))
# end

# function _get_or_load_cache(ds::GenomeDataset, idx::Int)
#     cached = _find_in_cache(ds._cache, idx)
#     !isnothing(cached) && return cached

#     genome = process_genome_one_side(ds.genome_dirs[idx]; side=ds.cds_side, pad=ds.pad)
#     _insert_cache!(ds, idx, genome)
#     return genome
# end

function Base.length(ds::GenomeDataset)
    return last(ds._cum_counts)
end

function Base.getindex(ds::GenomeDataset, idx::Int)::Tuple{Chromosome_X, Labels_y}
    idx in 1:length(ds) || throw(BoundsError(ds, idx))
    
    genome_idx = searchsortedfirst(ds._cum_counts, idx)
    fromcache = ds._cache[genome_idx]
    
    genome = if isnothing(fromcache)
        @warn "Accesing uncached genome $genome_idx"
        process_genome_one_side(ds.genome_dirs[genome_idx]; side=ds.cds_side, pad=ds.pad)
    else 
        fromcache
    end

    prev_total = genome_idx == 1 ? 0 : ds._cum_counts[genome_idx - 1]
    chrom_idx = idx - prev_total
    
    return (genome[1][chrom_idx], genome[2][chrom_idx])
end

function Base.lastindex(ds::GenomeDataset)
    return length(ds)
end

function Base.iterate(ds::GenomeDataset, state::Int)
    if state > length(ds)
        return nothing
    else
        return (ds[state], state + 1)
    end
end

function Base.iterate(ds::GenomeDataset)
    return iterate(ds, 1)
end