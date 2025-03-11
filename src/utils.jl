using GFF3
using DataFrames
using PrettyTables
using PlotlyJS

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
    @info "All dependencies found"
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
    precision(CM::ConfusionMTR)::Float

Возвращает precision (точность) из ConfusionMTR.
"""
function precision(CM::ConfusionMTR)
    return CM.TP / (CM.TP + CM.FP)
end

"""
    recall(CM::ConfusionMTR)::Float64

Возвращает recall (полноту) из ConfusionMTR.
"""
function recall(CM::ConfusionMTR)
    return CM.TP / (CM.TP + CM.FN)
end

"""
    f1_score(CM::ConfusionMTR)::Float64

Вычисляет F1-score на основе precision и recall.
"""
function f1_score(CM::ConfusionMTR)
    prec = precision(CM)
    rec = recall(CM)
    return 2 * (prec * rec) / (prec + rec)
end

"""
    false_dr(CM::ConfusionMTR)::Float64

Возвращает False Discovery Rate (FDR) из ConfusionMTR.
"""
function false_dr(CM::ConfusionMTR)
    return CM.FP / (CM.TP + CM.FP)
end



"""
    open_gff(name::String)::Vector{GFF3.Record}

Читает GFF-файл и возвращает список записей.
Аргументы:
- `name`: Имя файла.
"""
open_gff(name::String) = open(GFF3.Reader, name) do gff collect(gff) end

"""
    locus(record::GFF3.Record)::String

Формирует строку локуса для GFF3 записи.
Пример: "chr1:100..200".
"""
locus(record::GFF3.Record) = "$(GFF3.seqid(record)):$(GFF3.seqstart(record))..$(GFF3.seqend(record))"

"""
    locus(seqid::String)::Function

Создает функцию для формирования локуса с заданным seqid.
Пример: `locus("chr1")(1:100)` → "chr1:1..100".
"""
locus(seqid::String) = interval::UnitRange{Int64} -> "$seqid:$(interval.start)..$(interval.stop)"

"""
    ranges_from_GFF_records(list::Vector{GFF3.Record})::Vector{UnitRange{Int}}

Извлекает диапазоны из списка GFF3 записей.
"""
ranges_from_GFF_records(list::Vector{GFF3.Record}) = [
    GFF3.seqstart(x):GFF3.seqend(x) for x in list
]

"""
    filter_gff_region(sequence_header::String; regiontype::String, strand::String="+", phase::Int=-1)::Function

Создает фильтр для GFF3 записей по параметрам.
Аргументы:
- `sequence_header`: Имя хромосомы.
- `regiontype`: Тип области (например, "gene").
- `strand`: Направление ("+", "-", или не задано).
- `phase`: Фаза (или -1 для игнорирования).
Возвращает:
- Функцию фильтрации.
"""
filter_gff_region(sequence_header::String; regiontype::String, strand::String="+", phase::Int=-1) = list::Vector{GFF3.Record} -> filter(
    x -> all([
        GFF3.featuretype(x) == regiontype,
        GFF3.seqid(x) == sequence_header,
        phase == -1 || GFF3.hasphase(x) && GFF3.phase(x) == phase,
        string(GFF3.strand(x)) == strand
    ]), list
)


"""
    feature_strand_table(gff_entries::Vector{GFF3.Record}, chrom::String)::DataFrame

Создает таблицу распределения фичей по направлениям на хромосоме.
Столбцы:
- `Feature`: Тип фичи.
- `Positive_Strand`: Количество на позитивной цепи.
- `Negative_Strand`: Количество на негативной цепи.
"""
function feature_strand_table(gff_entries::Vector{GFF3.Record}, chrom::String)::DataFrame
    features = gff_entries .|> GFF3.featuretype |> Set

    df = DataFrame(
        Feature = String[],
        Positive_Strand = Int[],
        Negative_Strand = Int[]
    )
    
    for ft in features
        pos_count = gff_entries |> filter_gff_region(chrom; strand="+", regiontype=ft) |> length
        neg_count = gff_entries |> filter_gff_region(chrom; strand="-", regiontype=ft) |> length
        push!(df, (ft, pos_count, neg_count))
    end
    sort!(df, [:Positive_Strand, :Negative_Strand], rev=true)
    return df
end

"""
    rand_locus(gff_entries::Vector{GFF3.Record}, chrom::String, feature::String; strand::String="+")::String

Выбирает случайный локус для заданной фичи и хромосомы.
Аргументы:
- `gff_entries`: Список записей GFF.
- `chrom`: Имя хромосомы.
- `feature`: Тип фичи.
- `strand`: Направление (по умолчанию "+").
Возвращает:
- Случайный локус в формате строки.
"""
function rand_locus(gff_entries::Vector{GFF3.Record}, chrom::String, feature::String; strand::String="+")::String
    region = gff_entries |> filter_gff_region(chrom; strand=strand, regiontype=feature) |> x->rand(x)
    return locus(region)
end

"""
    extract_regions(gff_entries::Vector{GFF3.Record}, chrom::String, regtype="gene")::Tuple{Vector{UnitRange{Int}}, Vector{UnitRange{Int}}}

Извлекает диапазоны для позитивной и негативной цепей заданного типа.
Аргументы:
- `gff_entries`: Список записей GFF.
- `chrom`: Имя хромосомы.
- `regtype`: Тип области (по умолчанию "gene").
Возвращает:
- Кортеж с диапазонами для позитивной и негативной цепей.
"""
function extract_regions(gff_entries::Vector{GFF3.Record}, chrom::String, regtype="gene")
    features = gff_entries .|> GFF3.featuretype |> Set
    regtype ∉ features && return UnitRange{Int}[]

    gff_pos = gff_entries |> filter_gff_region(chrom; strand="+", regiontype=regtype) |> ranges_from_GFF_records# |> merge_ranges
    gff_neg = gff_entries |> filter_gff_region(chrom; strand="-", regiontype=regtype) |> ranges_from_GFF_records# |> merge_ranges

    return gff_pos, gff_neg
end


"""
    merge_ranges(r1::Vector{UnitRange{Int}}, r2::Vector{UnitRange{Int}})::Vector{UnitRange{Int}}

Объединяет пересекающиеся диапазоны из двух списков.
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
    merge_ranges(vr::Vector{Vector{UnitRange{Int}}})::Vector{UnitRange{Int}}

Объединяет все диапазоны из списка списков.
"""
function merge_ranges(vr::Vector{T})::T where T <: Vector{UnitRange{Int64}}
    reduce(merge_ranges, vr)
end

"""
    merge_ranges(r::Vector{UnitRange{Int}})::Vector{UnitRange{Int}}

Объединяет пересекающиеся диапазоны в списке.
"""
merge_ranges(r::Vector{UnitRange{Int64}})::Vector{UnitRange{Int64}} = merge_ranges(r, r)

"""
    range_intersection_matrix(ranges::Vector{UnitRange{Int}})::Matrix{Int}

Создает матрицу пересечений диапазонов (1 если пересекаются, иначе 0).
"""
function range_intersection_matrix(ranges::Vector{UnitRange{Int}})
    isintersecting(r1, r2) = !(r1.stop < r2.start || r2.stop < r1.start)
    
    n = length(ranges)
    intersection_matrix = [Int(isintersecting(ranges[i], ranges[j])) for i in 1:n, j in 1:n]

    return intersection_matrix
end

"""
    intersections_per_gene(locuses::Vector{UnitRange{Int}})::Vector{Int}

Считает количество пересечений для каждой гены в списке.
"""
function intersections_per_gene(locuses::Vector{UnitRange{Int}})
    IM = range_intersection_matrix(locuses)
    intersections = vec(sum(IM, dims=1))
    return intersections
end

"""
    len_intersection_matrix(ranges::Vector{UnitRange{Int}})::Matrix{Int}

Создает матрицу длин пересечений между диапазонами.
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
    intersections_length_per_gene(locuses::Vector{UnitRange{Int}})::Vector{Int}

Считает суммарную длину пересечений для каждой гена.
"""
function intersections_length_per_gene(locuses::Vector{UnitRange{Int}})
    LIM = len_intersection_matrix(locuses)
    intersection_lengths = vec(sum(LIM, dims=1))
    return intersection_lengths
end


"""
    longest_subcontig(subinters::Vector{UnitRange{Int}}, inters::Vector{UnitRange{Int}})::Vector{UnitRange{Int}}

Находит наибольшие подинтервалы, полностью вложенные в заданные интервалы.
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
    intersection_lengths(ref::Vector{UnitRange{Int}}, item::Vector{UnitRange{Int}})::Vector{Int}

Считает длины пересечений между всеми парами диапазонов.
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
    intersection_ratios(ref::Vector{UnitRange{Int}}, item::Vector{UnitRange{Int}})::Vector{Float64}

Считает доли пересечений относительно длины справочных диапазонов.
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