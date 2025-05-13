include("../src/utils.jl")
include("../src/plots.jl")

using Plots: plot, plot!
using Statistics
using Trapz
using DelimitedFiles

function threshold_metrics(predictions::Vector{Float32}, ground_truth::Vector{Bool})
    thresholds = Iterators.flatten([logrange(1e-8, 0.9e-1, 50),  range(1e-1, 1, 50)]) |> collect
    n = length(predictions)
    
    pos_cnt = count(ground_truth)
    neg_cnt = n - pos_cnt
    
    accuracies = Float64[]
    false_positives = Int[]
    false_negatives = Int[]
    tprs = Float64[]
    tnrs = Float64[]
    precisions = Float64[]
    f1_scores = Float64[]
    
    
    for t in thresholds
        preds_bool = predictions .>= t
        
        correct = sum(preds_bool .== ground_truth)
        push!(accuracies, correct / n)
        
        TP = sum(preds_bool .& ground_truth)
        FP = sum(preds_bool .& .!ground_truth)
        FN = sum((.!preds_bool) .& ground_truth)
        TN = sum((.!preds_bool) .& .!ground_truth)
        
        push!(false_positives, FP)
        push!(false_negatives, FN)
        
        push!(tprs, pos_cnt > 0 ? TP/pos_cnt : 0)
        push!(tnrs, neg_cnt > 0 ? TN/neg_cnt : 0)
        
        precision = TP / (TP + FP + 1e-9)
        push!(precisions, precision)
        f1 = 2 * precision * tprs[end] / (precision + tprs[end] + 1e-9)
        push!(f1_scores, f1)
    end
    
    fpr = [fp / neg_cnt for fp in false_positives]

    sorted_indices = sortperm(thresholds, rev=true)
    sorted_fpr = fpr[sorted_indices]
    sorted_tpr = tprs[sorted_indices]
    sorted_precision = precisions[sorted_indices]

    auc = trapz(sorted_fpr, sorted_tpr)
    precision_recall_auc = trapz(sorted_tpr, sorted_precision)

    
    return Dict(
        "thresholds" => collect(thresholds),
        "accuracy" => accuracies,
        "false_positives" => false_positives,
        "false_negatives" => false_negatives,
        "tpr" => tprs,
        "tnr" => tnrs,
        "precision" => precisions,
        "f1_score" => f1_scores,
        "auc_roc" => auc,
        "fpr" => fpr,
        "auc_pr" => precision_recall_auc,
        "pos_cnt" => pos_cnt,
        "neg_cnt" => neg_cnt,
    )
end

function plot_threshold_metrics(predictions::Vector{Float32}, ground_truth::Vector{Bool})
    metrics = threshold_metrics(predictions, ground_truth)
    thresholds = metrics["thresholds"]
    tpr = metrics["tpr"]
    precision = metrics["precision"]
    f1_scores = metrics["f1_score"]
    auc_roc = metrics["auc_roc"]
    auc_pr = metrics["auc_pr"]

    pos_cnt = metrics["pos_cnt"]
    neg_cnt = metrics["neg_cnt"]
    pos_rate = pos_cnt / (pos_cnt + neg_cnt)

    plt1 = plot(thresholds, tpr, xlabel="Threshold", ylabel="Recall (TPR)", labelfontsize=8,
                title="Recall (TPR) vs Threshold", legend=false, titlefontsize=10,)

    plt2 = plot(thresholds, precision, xlabel="Threshold", ylabel="Precision", labelfontsize=8,
                title="Precision vs Threshold", legend=false, titlefontsize=10,)

    plt3 = plot(thresholds, f1_scores, xlabel="Threshold", ylabel="F1 Score", labelfontsize=8,
                title="F1 Score vs Threshold", legend=false, titlefontsize=10,)

    plt4 = plot(
        metrics["fpr"], tpr, 
        label="ROC Curve", 
        xlabel="FPR (False Positive Rate)", 
        ylabel="TPR (Recall)", 
        title="ROC Curve (AUC = $(round(auc_roc, digits=3)))", 
        titlefontsize=10, 
        labelfontsize=8,
        legend=:bottomright,
    )
    
    plot!(plt4, [0, 1], [0, 1], 
            line=(:dash, :gray), 
            label="Random (AUC = 0.5)")

    plt5 = plot(thresholds, metrics["tnr"], 
            xlabel="Threshold", ylabel="Specificity (TNR)",
            title="Specificity (TNR) vs Threshold", legend=false,
            titlefontsize=10, labelfontsize=8,)
    
    plt6 = plot(precision, tpr, zcolor=f1_scores,
            xlabel="Precision", ylabel="Recall (TPR)",
            title="Precision-Recall Trade-off (AUC-PR = $(round(auc_pr, digits=3)))",
            legend=false, titlefontsize=10, labelfontsize=8,)

    plot(plt1, plt2, plt3, plt4, plt5, plt6,
         layout=(2, 3), size=(1200, 900), dpi=150)
end

function get_regions(bool_vector::Vector{Bool})
    regions = Tuple{Int, Int}[]
    n = length(bool_vector)
    start_idx = nothing

    for i in 1:n
        if bool_vector[i]
            if start_idx === nothing
                start_idx = i
            end
        else
            if start_idx !== nothing
                push!(regions, (start_idx, i - 1))
                start_idx = nothing
            end
        end
    end

    if start_idx !== nothing
        push!(regions, (start_idx, n))
    end

    return regions
end

function write_gff3(real::Vector{Bool}, predicted::Vector{Bool}, filename::String; seqid="chromosome1")
    real_regions = get_regions(real)
    predicted_regions = get_regions(predicted)

    real_set = Set(real_regions)
    correct_regions = []
    false_positive_regions = []

    for pr in predicted_regions
        if pr in real_set
            push!(correct_regions, pr)
        else
            push!(false_positive_regions, pr)
        end
    end

    open(filename, "w") do io
        println(io, "##gff-version 3")

        # Запись реальных регионов
        for (i, (start, stop)) in enumerate(real_regions)
            line = "$seqid\tReal\tCDS\t$start\t$stop\t.\t+\t.\tID=cds_real_$i;Name=Real_CDS_$i"
            println(io, line)
        end

        # Запись корректных предсказаний
        for (i, (start, stop)) in enumerate(correct_regions)
            line = "$seqid\tcorrect\tCDS\t$start\t$stop\t.\t+\t.\tID=cds_correct_$i;Name=Correct_Pred_$i"
            println(io, line)
        end

        # Запись ложных срабатываний
        for (i, (start, stop)) in enumerate(false_positive_regions)
            line = "$seqid\tfalse_positive\tCDS\t$start\t$stop\t.\t+\t.\tID=cds_fp_$i;Name=False_Positive_$i"
            println(io, line)
        end
    end
end

device!(3)

modelfile = "DATA/saved_models/2025-05-09T21:12:27.852slower_longer_starts/CuDevice(1)_epoch_012.flux"
m = deserialize(modelfile);
gpu_model = gpu(m.model)

pseudomonadota_accs = "DATA/genomes/genomes/" .* readdlm("DATA/subsets/pseudomonadota.tsv", '\t')[:, 1] .|> string
random_bacillota = "DATA/genomes/genomes/GCF_000091405.1"

genome_dir = random_bacillota # pseudomonadota_accs[end-200]
isdir(genome_dir)
ds_pseudo_1 = GenomeDataset([genome_dir], cds_side=:starts, pad=m.pad)

chrom, labels_pos, labels_neg = first(ds_pseudo_1)

X = _transform_X(Flux.onehotbatch(chrom, ('A','C','G','T'), 'C')) |> gpu
y = labels_pos

ŷ_frommodel = @time cpu(gpu_model(X))

tm = threshold_metrics(ŷ_frommodel, y)
optimal_threshold = tm["thresholds"][argmax(tm["f1_score"])]
ŷ = ŷ_frommodel .> optimal_threshold

total_ends = sum(y)
predicted_ends = sum(ŷ)
matches_ends = sum(y .* ŷ)

best_f1 = maximum(tm["f1_score"])
plot_threshold_metrics(ŷ_frommodel, y)
plot_mismatch_histograms(m.cds_ranges_true, m.cds_ranges_model, m.fp_shifts)

sequence = join(chrom[m.pad+1:end-m.pad])

write_gff3(y, Vector{Bool}(ŷ), "chromosome_bacillota.gff3", seqid="chr1")
FastaIO.writefasta("chromosome_bacillota.fa", [("chr1", sequence)])
