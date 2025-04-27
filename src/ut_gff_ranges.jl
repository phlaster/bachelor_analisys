using GFF3
using DataFrames

"""
    open_gff(name)

Read and parse a GFF3 file, returning all records as a vector.
"""
open_gff(name::T) where T <: AbstractString = open(GFF3.Reader, name) do gff collect(gff) end

"""
    locus(record)

Create a standardized locus string from a GFF3 record in "seqid:start..end" format. 
Example: "chrX:5000..7500" represents a feature spanning positions 5000 to 7500 on chromosome X.
"""
locus(record::GFF3.Record) = "$(GFF3.seqid(record)):$(GFF3.seqstart(record))..$(GFF3.seqend(record))"

"""
    locus(seqid)

Create a locus string generator for a specific chromosome/sequence identifier. The returned 
function accepts a UnitRange{Int} to produce locus strings (e.g., locus("chr2")(1000:2000) → "chr2:1000..2000").
"""
locus(seqid::T) where T <: AbstractString = interval::UnitRange{Int64} -> "$seqid:$(interval.start)..$(interval.stop)"

"""
    ranges_from_GFF_records(list[, side={:left|:right})

Extract genomic position ranges (or borders of a chosen side) from GFF3 records.

```
julia> ranges_from_GFF_records(some_gff_list)
2-element Vector{UnitRange{Int64}}:
 1:15
 22:30

julia> ranges_from_GFF_records(some_gff_list, :left)
2-element Vector{Int64}:
  1
 22

julia> ranges_from_GFF_records(some_gff_list, :right)
2-element Vector{Int64}:
 15
 30
````
"""
ranges_from_GFF_records(list::Vector{GFF3.Record}) = [
    GFF3.seqstart(x):GFF3.seqend(x) for x in list
]

function ranges_from_GFF_records(list::Vector{GFF3.Record}, side::Symbol)
    if side == :left
        return [GFF3.seqstart(x) for x in list]
    elseif side == :right
        return [GFF3.seqend(x) for x in list]
    else
        return Int[]
    end
end

function ranges_from_GFF_records!(chromosome::Vector{Bool}, list::Vector{GFF3.Record}, side::Symbol)
    isempty(list) || isempty(chromosome) && return
    L = length(chromosome)
    if side == :left
        @simd for x in list
            position = mod1(GFF3.seqstart(x), L)
            chromosome[position] |= true
        end
    elseif side == :right
        @simd for x in list
            position = mod1(GFF3.seqend(x), L)
            chromosome[position] |= true
        end
    end
end


"""
    filter_gff_region(; sequence_header, regiontype, strand, phase)

Create a filter function for GFF3 records based on multiple criteria:
- `sequence_header`: Chromosome/contig name to match, if empty choses all;
- `regiontype`: Feature type (e.g., "gene", "exon"), if empty choses all;
- `strand`: Strand direction ("+" or "-"), if empty choses both;
- `phase`: Translation phase (0-2), use -1 to disable phase filtering;
The returned function filters vectors of GFF3.Record objects.
"""
filter_gff_region(;
    sequence_header="",
    regiontype="",
    strand="",
    phase::Int=-1
) = list::Vector{GFF3.Record} -> filter(
    x -> all([
        isempty(regiontype)      || GFF3.featuretype(x) == regiontype,
        isempty(sequence_header) || GFF3.seqid(x) == sequence_header,
        phase == -1              || GFF3.hasphase(x) && GFF3.phase(x) == phase,
        isempty(strand)          || string(GFF3.strand(x)) == strand
    ]), list
)


"""
    feature_strand_table(gff_entries, chrom)

Generate a summary table showing strand distribution of genomic features. Returns a DataFrame 
with columns: Feature type, Positive Strand count, and Negative Strand count for features 
on the specified chromosome.
"""
function feature_strand_table(gff_entries::Vector{GFF3.Record}, chrom="")
    features = gff_entries .|> GFF3.featuretype |> Set

    df = DataFrame(
        Feature = String[],
        Positive_Strand = Int[],
        Negative_Strand = Int[]
    )
    
    for ft in features
        pos_count = gff_entries |> filter_gff_region(; sequence_header=chrom, strand="+", regiontype=ft) |> length
        neg_count = gff_entries |> filter_gff_region(; sequence_header=chrom, strand="-", regiontype=ft) |> length
        push!(df, (ft, pos_count, neg_count))
    end
    DataFrames.sort!(df, [:Positive_Strand, :Negative_Strand], rev=true)
    return df
end

"""
    rand_locus(gff_entries, chrom, feature; strand)

Randomly select a genomic locus from features of specified type and strand. Useful for 
sampling representative regions from annotation data.
"""
function rand_locus(gff_entries::Vector{GFF3.Record}, chrom::T, feature::T; strand::T="")::String where T <: AbstractString
    region = gff_entries |> filter_gff_region(; sequence_header=chrom, strand=strand, regiontype=feature) |> x->rand(x)
    return locus(region)
end

"""
    extract_regions(gff_entries, chrom, regtype)

Separate genomic features into positive and negative strand ranges. Returns a tuple 
containing two vectors of UnitRanges: (positive_strand_ranges, negative_strand_ranges) 
for the specified feature type on the given chromosome.
"""
function extract_regions(gff_entries::Vector{GFF3.Record}, chrom::T, regtype="gene") where T <: AbstractString
    features = gff_entries .|> GFF3.featuretype |> Set
    regtype ∉ features && return UnitRange{Int}[]

    gff_pos = gff_entries |> filter_gff_region(; sequence_header=chrom, strand="+", regiontype=regtype) |> ranges_from_GFF_records# |> merge_ranges
    gff_neg = gff_entries |> filter_gff_region(; sequence_header=chrom, strand="-", regiontype=regtype) |> ranges_from_GFF_records# |> merge_ranges

    return gff_pos, gff_neg
end


"""
    merge_ranges(r1, r2)

Combine and simplify overlapping or adjacent genomic ranges. Takes two vectors of ranges 
and returns a single vector of non-overlapping ranges that cover all original positions.
"""
function merge_ranges(r1::T, r2::T)::T where T <: Vector{UnitRange{Int64}}
    isempty(r1) && return r2
    isempty(r2) && return r1
    
    combined_ranges = vcat(r1, r2)
    length(combined_ranges) == 1 && return combined_ranges

    sorted_ranges = sort(combined_ranges, by=x -> x.start)
    merged_ranges = [sorted_ranges[1]]

    for current_range in sorted_ranges[2:end]
        last_merged_range = merged_ranges[end]

        if current_range.start <= last_merged_range.stop + 1
            new_range = last_merged_range.start:max(last_merged_range.stop, current_range.stop)
            merged_ranges[end] = new_range
        else
            push!(merged_ranges, current_range)
        end
    end

    return merged_ranges
end

"""
    merge_ranges(vr)

Combine multiple sets of genomic ranges into a single simplified set. Accepts a vector 
of range vectors and returns merged ranges through successive pairwise merging.
"""
function merge_ranges(vr::Vector{T})::T where T <: Vector{UnitRange{Int64}}
    reduce(merge_ranges, vr)
end

merge_ranges(r::Vector{UnitRange{Int64}})::Vector{UnitRange{Int64}} = merge_ranges(r, r)

"""
    range_intersection_matrix(ranges)

Create an adjacency matrix showing pairwise range overlaps. Matrix elements are 1 if ranges 
i and j overlap (or touch), 0 otherwise. Diagonal represents self-comparisons.
"""
function range_intersection_matrix(ranges::Vector{UnitRange{Int}})
    isintersecting(r1, r2) = !(r1.stop < r2.start || r2.stop < r1.start)
    
    n = length(ranges)
    intersection_matrix = [Int(isintersecting(ranges[i], ranges[j])) for i in 1:n, j in 1:n]

    return intersection_matrix
end

"""
    intersections_per_gene(locuses)

Calculate overlap frequency for genomic features. Returns a vector where each element 
indicates how many other ranges intersect with the corresponding input range.
"""
function intersections_per_gene(locuses::Vector{UnitRange{Int}})
    IM = range_intersection_matrix(locuses)
    intersections = vec(sum(IM, dims=1))
    return intersections
end

"""
    len_intersection_matrix(ranges)

Create a matrix of basepair overlap lengths between genomic ranges. Only calculates 
upper triangle values to avoid duplicate computations between range pairs.
"""
function len_intersection_matrix(ranges::Vector{UnitRange{Int}})
    n = length(ranges)
    
    intersection_matrix = zeros(Int, n, n)

    @inbounds for i in 1:n
        @simd for j in i+1:n
            intersection_matrix[i, j] = length(intersect(ranges[i], ranges[j]))
        end
    end

    return intersection_matrix
end

"""
    intersections_length_per_gene(locuses)

Calculate total overlap length per genomic feature. Returns a vector where each element 
represents the sum of all basepair overlaps between a range and others in the set.
"""
function intersections_length_per_gene(locuses::Vector{UnitRange{Int}})
    LIM = len_intersection_matrix(locuses)
    intersection_lengths = vec(sum(LIM, dims=1))
    return intersection_lengths
end

"""
    longest_subcontig(subinters, inters)

Find maximal continuous subregions within reference intervals. For each interval in `inters`, 
identifies the largest contiguous span covered by ranges from `subinters` that are fully 
contained within it.
"""
function longest_subcontig(subinters::Vector{UnitRange{Int}}, inters::Vector{UnitRange{Int}})
    result = Vector{UnitRange{Int}}()
    
    for range2 in inters
        sub_ranges = filter(range1 -> range1.start >= range2.start && range1.stop <= range2.stop, subinters)
        if !isempty(sub_ranges)
            min_start = minimum(range1.start for range1 in sub_ranges)
            max_stop = maximum(range1.stop for range1 in sub_ranges)            
            push!(result, min_start:max_stop)
        end
    end
    
    return result
end

"""
    intersection_lengths(ref, item)

Calculate pairwise overlap lengths between two sets of genomic ranges. Returns all 
intersection lengths between ranges from `ref` and `item`, including multiple overlaps 
for nested features.
"""
function intersection_lengths(ref::Vector{UnitRange{Int}}, item::Vector{UnitRange{Int}})
    result = Int[]
    
    for r1 in ref
        for r2 in item
            intersection_start = max(r1.start, r2.start)
            intersection_stop = min(r1.stop, r2.stop)
            if intersection_start <= intersection_stop
                push!(result, intersection_stop - intersection_start + 1)
            end
        end
    end
    
    return result
end


"""
    intersection_ratios(ref, item)

Calculate proportional overlaps between reference and query ranges. Returns ratios of 
intersection lengths relative to reference range lengths, useful for measuring coverage 
efficiency.
"""
function intersection_ratios(ref::Vector{UnitRange{Int}}, item::Vector{UnitRange{Int}})
    result = Float64[]
    
    for r1 in ref
        for r2 in item
            intersection_start = max(r1.start, r2.start)
            intersection_stop = min(r1.stop, r2.stop)
            if intersection_start <= intersection_stop
                reflen = length(r1)
                ilen = intersection_stop - intersection_start + 1
                push!(result, ilen/reflen)
            end
        end
    end
    
    return result
end