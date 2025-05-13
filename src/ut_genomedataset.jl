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
        @assert cds_side in [:starts, :stops, :inner] "wrong side symbol :$side, must be either :starts or :stops or :inner"
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

function CDS_borders_both_strands(gff_data, chrom_names, chrom_lengths, side::Symbol)
    @assert !isempty(gff_data) "Empty gff data"
    @assert !isempty(chrom_names) "Empty chromosomes list"
    
    preallocated_pos = [
        zeros(Bool, chrom_length)
        for chrom_length in chrom_lengths
    ]
    preallocated_neg = [
        zeros(Bool, chrom_length)
        for chrom_length in chrom_lengths
    ]
    for (chrom_name, labels_pos, labels_neg) in zip(chrom_names, preallocated_pos, preallocated_neg)
        cds_filter_pos = filter_gff_region(;
            sequence_header=chrom_name,
            regiontype="CDS",
            strand="+"
        )
        cds_filter_neg = filter_gff_region(;
            sequence_header=chrom_name,
            regiontype="CDS",
            strand="-"
        )

        CDS_regions_pos = cds_filter_pos(gff_data)
        CDS_regions_neg = cds_filter_neg(gff_data)

        if side == :starts
            ranges_from_GFF_records!(labels_pos, CDS_regions_pos, :left)
            ranges_from_GFF_records!(labels_neg, CDS_regions_neg, :right)
        elseif side == :stops
            ranges_from_GFF_records!(labels_pos, CDS_regions_pos, :right)
            ranges_from_GFF_records!(labels_neg, CDS_regions_neg, :left)
        elseif side == :inner
            ranges_from_GFF_records!(labels_pos, CDS_regions_pos, :inner)
            ranges_from_GFF_records!(labels_neg, CDS_regions_neg, :inner)
        else
            throw(ArgumentError)
        end
    end
    return preallocated_pos, preallocated_neg
end

function process_genome_one_side(genome_dir::T; side::Symbol=:starts, pad::Int=0) where T <: AbstractString
    @assert side in [:starts, :stops, :inner] "wrong side symbol: $side, must be either :starts or :stops or :inner"
    @assert isdir(genome_dir) "No such genomic directory"
    file_prefix = last(splitpath(genome_dir))
    fastaname = joinpath(genome_dir, "$file_prefix.fna")
    gffname = joinpath(genome_dir, "$file_prefix.gff")
    @assert all(isfile.([fastaname, gffname])) "Genomic files not found"

    
    fasta_data = readfasta(fastaname)
    fasta_sequences = getindex.(fasta_data, 2)
    dna_padded = [add_pad(uppercase.(collect(fs)), pad) for fs in fasta_sequences]
    # dna_encoded =  Flux.onehotbatch.(dna_padded, Ref(('A', 'C', 'G', 'T')), 'C')

    
    chrom_names = getindex.(fasta_data, 1) .|> split .|> first
    chrom_lengths = length.(fasta_sequences)
    gff_data = open_gff(gffname)
    
    labels_pos, labels_neg = CDS_borders_both_strands(gff_data, chrom_names, chrom_lengths, side)
    
    return (dna_padded, labels_pos, labels_neg)
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