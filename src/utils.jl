using GFF3
using DataFrames
using PrettyTables
using FastaIO

function clear_last_lines(n::Int)
    for _ in 1:n
        print("\e[1A")
        print("\e[2K")
    end
    flush(stdout)
end

function check_dependencies(execs)
    for cmd in execs
        try
            run(pipeline(`which $cmd`, devnull))
        catch
            @error "Required command '$cmd' not found in PATH"
            exit(1)
        end
    end
end

function wait_tasks(tasks::Vector{Task})
    while true
        for t in tasks
            if istaskdone(t)
                return t
            end
        end
        sleep(0.5)
    end
end

"""
    ConfusionMTR

Структура для представления матрицы результативности и связанных метрик классификации.

Поля:
- `TP::Int`: True Positives
- `FP::Int`: False Positives
- `TN::Union{Int, Nothing}`: True Negatives (может быть `nothing` для некоторых методов)
- `FN::Int`: False Negatives
- `prec::Float64`: Precision
- `rec::Float64`: Recall
- `f1::Float64`: F1-score
- `fdr::Float64`: False Discovery Rate
- `mtype::String`: Тип матрицы (например, "Confusion" или "Contingency")
"""
struct ConfusionMTR
    TP::Int
    FP::Int
    TN::Union{Int, Nothing}
    FN::Int

    prec::Float64
    rec::Float64
    f1::Float64
    fdr::Float64

    mtype::String

    """
        ConfusionMTR(mtype::String, cm::NamedTuple)

    Конструктор для инициализации структуры на основе именованного кортежа с TP, FP, TN, FN.
    Вычисляет метрики precision, recall, F1 и FDR.
    """
    function ConfusionMTR(mtype::String, cm::@NamedTuple{TP::Int64, FP::Int64, TN::Int64, FN::Int64})
        prec = cm.TP / (cm.TP + cm.FP)
        rec = cm.TP / (cm.TP + cm.FN)
        f1 = 2 * (prec * rec) / (prec + rec)
        fdr = cm.FP / (cm.TP + cm.FP)
        new(cm.TP, cm.FP, cm.TN, cm.FN, prec, rec, f1, fdr, mtype)
    end
end

"""
    ConfusionMTR(ref_r::Vector{UnitRange{Int}}, test_r::Vector{UnitRange{Int}})

Создает ConfusionMTR для сравнения двух наборов позиций в геноме.
Аргументы:
- `ref_r`: Справочные диапазоны.
- `test_r`: Тестовые диапазоны.
Возвращает:
- ConfusionMTR с вычисленными TP, FP, TN, FN и метриками.
"""
function ConfusionMTR(ref_r::T, test_r::T) where T <: Vector{UnitRange{Int64}}
    reference_positions = Set{Int}()
    for range in ref_r
        union!(reference_positions, range)
    end

    test_positions = Set{Int}()
    for range in test_r
        union!(test_positions, range)
    end

    TP = length(reference_positions ∩ test_positions)
    FP = length(setdiff(test_positions, reference_positions))
    FN = length(setdiff(reference_positions, test_positions))

    all_positions = Set{Int}(collect(1:maximum(reference_positions ∪ test_positions; init=0)))

    TN = length(setdiff(all_positions, reference_positions ∪ test_positions))

    return ConfusionMTR("Confusion", (TP=TP, FP=FP, TN=TN, FN=FN))
end

"""
    contingency(ref::Vector{UnitRange{Int}}, predicted::Vector{UnitRange{Int}}; threshold=0.5) -> ConfusionMTR

Вычисляет Contingency матрицу с учетом порога перекрытия.
Аргументы:
- `ref`: Справочные геновые диапазоны.
- `predicted`: Предсказанные диапазоны.
- `threshold`: Минимальная доля перекрытия (по умолчанию 0.5).
Возвращает:
- ConfusionMTR с TP, FP, FN (TN всегда -1).
"""
function contingency(ref::Vector{UnitRange{Int}}, predicted::Vector{UnitRange{Int}}; threshold=0.5)
    TP = 0
    FP = 0
    FN = 0
    
    for pred in predicted
        found_match = false
        
        for ref_gene in ref
            overlap_length = length(intersect(pred, ref_gene))
            min_length = min(length(pred), length(ref_gene))
            
            if overlap_length / min_length >= threshold
                TP += 1
                found_match = true
                break
            end
        end
        
        if !found_match
            FP += 1
        end
    end
    
    for ref_gene in ref
        found_match = false
        
        for pred in predicted
            overlap_length = length(intersect(pred, ref_gene))
            min_length = min(length(pred), length(ref_gene))
            
            if overlap_length / min_length >= threshold
                found_match = true
                break
            end
        end
        
        if !found_match
            FN += 1
        end
    end
    
    return ConfusionMTR("Contigency", (TP=TP, FP=FP, TN=-1, FN=FN))
end


function Base.show(io::IO, cm::ConfusionMTR)
    println(io, "$(cm.mtype) Matrix:")
    data = [
        "Real +" cm.TP        cm.FN                     cm.TP+cm.FN
        "Real -" cm.FP        ifelse(cm.TN>0,cm.TN,"-") cm.FP+max(0,cm.TN)
        "Total"  cm.TP+cm.FP  cm.FN+max(0,cm.TN)        cm.TP+cm.FP+cm.FN+max(0,cm.TN)
        ]
    pt = pretty_table(String, data; header = [" ", "Pred +", "Pred -", "Total"])
    println(io, pt)
    println(io, "Prescision = $(round(cm.prec, digits=3))")
    println(io, "Recall     = $(round(cm.rec, digits=3))")
    println(io, "F1         = $(round(cm.f1, digits=3))")
    println(io, "FDR        = $(round(cm.fdr, digits=3))")
end

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


"""
    filter_gff_region(; sequence_header, regiontype, strand, phase)

Create a filter function for GFF3 records based on multiple criteria:
- `sequence_header`: Chromosome/contig name to match, if empty choses all
- `regiontype`: Required feature type (e.g., "gene", "exon"), 
- `strand`: Strand direction ("+" or "-"), if empty choses all
- `phase`: Translation phase (0-2), use -1 to disable phase filtering
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
function feature_strand_table(gff_entries::Vector{GFF3.Record}, chrom::T) where T <: AbstractString
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
    sort!(df, [:Positive_Strand, :Negative_Strand], rev=true)
    return df
end

"""
    rand_locus(gff_entries, chrom, feature; strand)

Randomly select a genomic locus from features of specified type and strand. Useful for 
sampling representative regions from annotation data.
"""
function rand_locus(gff_entries::Vector{GFF3.Record}, chrom::T, feature::T; strand::T="+")::String where T <: AbstractString
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


function onehot_encode_dna(seq::AbstractString; revcomp=false)
    N = length(seq)
    encoding = zeros(UInt8, N, 4)

    if revcomp
        comp_mapping = Dict{Char,UInt8}('A'=>0x1, 'C'=>0x2, 'G'=>0x3, 'T'=>0x4)
        for (i, ch) in enumerate(Iterators.reverse(uppercase(seq)))
            if haskey(comp_mapping, ch)
                col = comp_mapping[ch]
                encoding[i, col] = 0x1
            end
        end
    else
        mapping = Dict{Char,UInt8}('T'=>0x1, 'G'=>0x2, 'C'=>0x3, 'A'=>0x4)
        for (i, ch) in enumerate(uppercase(seq))
            if haskey(mapping, ch)
                col = mapping[ch]
                encoding[i, col] = 0x1
            end
        end
    end
    return encoding
end


function CDS_borders_in_GFF(gff_data, chromosomes::Vector{<:AbstractString}; strand="+")
    @assert !isempty(gff_data) "Empty gff data"
    @assert !isempty(chromosomes) "Empty chromosomes list"
    @assert strand in ["+", "-"] "wrong strand identifier: $strand"
    
    cds_starts = Vector{UInt8}[]
    cds_stops = Vector{UInt8}[]
    for chrom in chromosomes
        region_filter = filter_gff_region(; sequence_header=chrom, regiontype="region")

        chrom_length = gff_data |> region_filter |> ranges_from_GFF_records |> only |> maximum

        cds_filter = filter_gff_region(;
            sequence_header=chrom,
            regiontype="CDS",
            strand=strand
        )
        intervals = ranges_from_GFF_records(cds_filter(gff_data))
        

        chrom_starts = zeros(UInt8, chrom_length)
        chrom_stops = zeros(UInt8, chrom_length)

        if isempty(intervals)
            push!(cds_starts, chrom_starts)
            push!(cds_stops, chrom_stops)
            continue
        end

        starts = replace(
            getfield.(intervals, :start) .% chrom_length,
            0 => chrom_length
        )
        stops = replace(
            getfield.(intervals, :stop) .% chrom_length,
            0 => chrom_length
        )
        @assert length(starts) == length(stops) "start-stop count mismatch"
        chrom_starts[starts] .= 0x1
        chrom_stops[stops] .= 0x1
        
        push!(cds_starts, chrom_starts)
        push!(cds_stops, chrom_stops)
    end
    
    return strand == "+" ?
        (cds_starts, cds_stops) :
        (reverse!.(cds_starts), reverse!.(cds_stops))
end

function digitize_genome(genome_dir::T) where T <: AbstractString
    @assert isdir(genome_dir) "No such genomic directory"
    file_prefix = last(splitpath(genome_dir))
    fastaname = joinpath(genome_dir, "$file_prefix.fna")
    gffname = joinpath(genome_dir, "$file_prefix.gff")
    @assert all(isfile.([fastaname, gffname])) "Genomic files not found"
    
    fasta_data = readfasta(fastaname)
    fasta_strings = getindex.(fasta_data, 2)
    
    dna_matrix_pos = onehot_encode_dna.(fasta_strings; revcomp=false)
    dna_matrix_neg = onehot_encode_dna.(fasta_strings; revcomp=true)
    
    chromosomes = getindex.(fasta_data, 1) .|> split .|> first
    gff_data = open_gff(gffname)
    
    starts_pos, ends_pos = CDS_borders_in_GFF(gff_data, chromosomes; strand="+")
    starts_neg, ends_neg = CDS_borders_in_GFF(gff_data, chromosomes; strand="-")

    @assert all(==(length(fasta_strings)), [
        length(dna_matrix_pos),
        length(dna_matrix_neg),
        length(starts_pos),
        length(ends_pos),
        length(starts_neg),
        length(ends_neg)
    ]) "Dimentions do not match"

    return dna_matrix_pos, dna_matrix_neg, starts_pos, ends_pos, starts_neg, ends_neg
end