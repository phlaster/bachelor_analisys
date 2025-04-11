using PrettyTables


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
    specificity::Float64
    support::Int


    mtype::String

    """
        ConfusionMTR(mtype::String, cm::NamedTuple)

    Конструктор для инициализации структуры на основе именованного кортежа с TP, FP, TN, FN.
    Вычисляет метрики precision, recall, F1 и FDR.
    """
    function ConfusionMTR(mtype::String, cm::@NamedTuple{TP::Int64, FP::Int64, TN::Int64, FN::Int64})
        prec = (cm.TP + cm.FP) == 0 ? 0.0 : cm.TP / (cm.TP + cm.FP)
        rec = (cm.TP + cm.FN) == 0 ? 0.0 : cm.TP / (cm.TP + cm.FN)
        f1 = 2 * (prec + rec) == 0 ? 0.0 : 2 * prec * rec / (prec + rec)
        fdr = (cm.TP + cm.FP) == 0 ? 0.0 : cm.FP / (cm.TP + cm.FP)
        specificity = (cm.TN + cm.FP) == 0 ? 0.0 : cm.TN / (cm.TN + cm.FP)
        support = cm.TP + cm.FN
        new(cm.TP, cm.FP, cm.TN, cm.FN, prec, rec, f1, fdr, specificity, support, mtype)
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
    println(io, "Specificity= $(round(cm.specificity, digits=3))")
    println(io, "Support    = $(cm.support)")
end

function multiclass_confusion(ground_truth, prediction;
    classes=sort!(unique(Iterators.flatten((ground_truth, prediction))))
)
    if length(ground_truth) != length(prediction)
        throw(DimensionMismatch("The input vectors must have the same length."))
    end

    n_classes = length(classes)    
    cm = zeros(Int, (n_classes, n_classes))
    class_to_idx = Dict{Int, Int}(class => idx for (idx, class) in enumerate(classes))
    @inbounds for (actual, predicted) in zip(ground_truth, prediction)
        row = class_to_idx[actual]
        col = class_to_idx[predicted]
        cm[row, col] += 1
    end
    return cm, classes
end

function metrics_from_cm(cm, classes)
    n_classes = length(classes)
    if size(cm) != (n_classes, n_classes)
        throw(DimensionMismatch("Confusion matrix must be a square matrix matching the number of classes."))
    end
    
    metrics_dict = Dict{Int, NamedTuple}()
    total_samples = sum(cm)
    
    for (i, class) in enumerate(classes)
        tp = cm[i, i]
        fp = sum(cm[:, i]) - tp
        fn = sum(cm[i, :]) - tp
        support = sum(cm[i, :])
        tn = total_samples - (tp + fp + fn)
        
        prec_denom = tp + fp
        prec = prec_denom == 0 ? NaN : tp / prec_denom
        
        rec_denom = tp + fn
        rec = rec_denom == 0 ? NaN : tp / rec_denom
        
        f1_denom = prec + rec
        f1 = f1_denom == 0 ? NaN : 2 * (prec * rec) / f1_denom
        
        fdr_denom = tp + fp
        fdr = fdr_denom == 0 ? NaN : fp / fdr_denom
        
        spec_denom = tn + fp
        specificity = spec_denom == 0 ? NaN : tn / spec_denom
        
        metrics = (
            prec = prec,
            rec = rec,
            f1 = f1,
            fdr = fdr,
            specificity = specificity,
            support = support
        )
        
        metrics_dict[class] = metrics
    end
    
    return metrics_dict, cm, classes
end