#!/usr/bin/env julia

const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)

using LuxCUDA
using Flux
using MLUtils
using Statistics
using ProgressMeter

# Структура загрузчика данных для одной хромосомы/стренда
struct ChromosomeDataLoader
    dna::Matrix{Float32}      # (4, seq_len) - one-hot encoding ДНК
    borders::Vector{Float32}  # (seq_len,) - метки начала генов
    window_size::Int          # Размер окна
    step::Int                 # Шаг между окнами
    positions::Vector{Int}    # Позиции для извлечения окон
end

# Конструктор загрузчика
function ChromosomeDataLoader(dna::Matrix{Float32}, borders::Vector{Float32}, window_size::Int, step::Int)
    seq_len = size(dna, 2)
    @assert size(dna, 1) == 4 "DNA matrix must have 4 rows"
    @assert length(borders) == seq_len "Borders length must match sequence length"
    positions = collect(1:step:seq_len)
    return ChromosomeDataLoader(dna, borders, window_size, step, positions)
end

# Извлечение окна с учетом кольцевости
function get_window(loader::ChromosomeDataLoader, pos::Int, buffer_dna::Matrix{Float32}, buffer_borders::Vector{Float32})
    seq_len = size(loader.dna, 2)
    @views for k in 1:loader.window_size
        idx = mod1(pos + k - 1, seq_len)
        buffer_dna[:, k] = loader.dna[:, idx]
        buffer_borders[k] = loader.borders[idx]
    end
    return buffer_dna, buffer_borders
end

# Создание всех окон из всех хромосом/стрендов
function create_all_windows(all_loaders::Vector{ChromosomeDataLoader})
    all_data = Vector{Tuple{Matrix{Float32}, Vector{Float32}}}()
    for loader in all_loaders
        buffer_dna = Matrix{Float32}(undef, 4, loader.window_size)
        buffer_borders = Vector{Float32}(undef, loader.window_size)
        for pos in loader.positions
            X, y = get_window(loader, pos, buffer_dna, buffer_borders)
            push!(all_data, (copy(X), copy(y)))
        end
    end
    return all_data
end

# Сборка батчей
function collate(batch)
    X = cat([b[1] for b in batch]..., dims=3)  # (4, window_size, batch_size)
    X = permutedims(X, (2, 1, 3))              # (window_size, 4, batch_size)
    y = cat([b[2] for b in batch]..., dims=2)  # (window_size, batch_size)
    return X, y
end

# Создание модели
function create_model()
    Chain(
        Conv((5,), 4 => 32, relu, pad=SamePad()),  # (window_size, 4, batch_size) -> (window_size, 32, batch_size)
        Conv((5,), 32 => 64, relu, pad=SamePad()), # -> (window_size, 64, batch_size)
        Conv((1,), 64 => 1, identity, pad=SamePad()), # -> (window_size, 1, batch_size)
        x -> sigmoid.(x),
        x -> dropdims(x, dims=2)  # -> (window_size, batch_size)
    ) |> gpu
end

# Функция потерь
loss(y_pred, y) = mean(Flux.binarycrossentropy.(y_pred, y))

function weighted_loss(y_pred, y, weight_pos=10.0)
    bce = Flux.binarycrossentropy.(y_pred, y)  # Бинарная кросс-энтропия
    weights = y .* weight_pos .+ (1 .- y)      # Вес weight_pos для y=1, 1 для y=0
    mean(bce .* weights)                       # Средняя взвешенная потеря
end

# Функция для вычисления метрик на тестовых данных
function evaluate_model(model, test_dataloader)
    total_tp = 0
    total_fp = 0
    total_fn = 0
    threshold = 0.5  # Порог для бинаризации предсказаний

    for (X, y) in test_dataloader
        X_gpu = X |> gpu
        y_gpu = y |> gpu
        y_pred = model(X_gpu) .> threshold  # Бинарные предсказания
        y_true = y_gpu .== 1

        tp = sum(y_pred .& y_true)  # True positives
        fp = sum(y_pred .& .!y_true)  # False positives
        fn = sum(.!y_pred .& y_true)  # False negatives

        total_tp += tp
        total_fp += fp
        total_fn += fn
    end

    precision = total_tp / (total_tp + total_fp + eps(Float32))  # Добавляем eps для избежания деления на 0
    recall = total_tp / (total_tp + total_fn + eps(Float32))
    f1 = 2 * precision * recall / (precision + recall + eps(Float32))

    return precision, recall, f1
end

# Обучение модели с оценкой на тестовых данных
function train_model!(model, state, train_dataloader, test_dataloader, epochs=10, weight_pos=10.0)
    for epoch in 1:epochs
        total_loss = 0.0
        num_batches = 0
        @showprogress desc="Progress for epoch $epoch" for (X, y) in train_dataloader
            X_gpu = X |> gpu
            y_gpu = y |> gpu
            gs = gradient(model -> weighted_loss(model(X_gpu), y_gpu, weight_pos), model)
            Flux.update!(state, model, gs[1])
            total_loss += weighted_loss(model(X_gpu), y_gpu, weight_pos)
            num_batches += 1
        end
        avg_loss = total_loss / num_batches
        println("Эпоха $epoch, средний loss: $avg_loss")

        # Оценка на тестовых данных
        precision, recall, f1 = evaluate_model(model, test_dataloader)
        println("Тестовые метрики: Precision=$precision, Recall=$recall, F1=$f1")
    end
end

# Основная функция
function main()
    train_directories = "DATA/genomes/genomes/" .* readdir("DATA/genomes/genomes")[1:10]
    test_directories = "DATA/genomes/genomes/" .* readdir("DATA/genomes/genomes")[11:15]
    window_size = 500
    step = 50
    batch_size = 32
    epochs = 3

    # Счётчики классов для взвешивания функции потерь
    pos_count = 0
    neg_count = 0

    # Создаем загрузчики для тренировочных данных
    train_loaders = ChromosomeDataLoader[]
    @showprogress desc="Loading train genomes" for dir in train_directories
        d = digitize_genome_one_side(dir)
        all_dna = vcat(d.dna_pos, d.dna_neg)
        all_borders = vcat(d.borders_pos, d.borders_neg)
        for (dna, borders) in zip(all_dna, all_borders)
            loader = ChromosomeDataLoader(dna, borders, window_size, step)
            push!(train_loaders, loader)
        end
        pos_count += sum(sum.(all_borders))  # Число положительных примеров
        neg_count += sum(length.(all_borders)) - pos_count  # Число отрицательных примеров
    end
    weight_pos = neg_count / pos_count  # Вес для положительного класса

    # Собираем все окна из тренировочных данных
    train_data = create_all_windows(train_loaders)
    train_dataloader = Flux.DataLoader(train_data, batchsize=batch_size, shuffle=true, collate=collate)

    # Создаем загрузчики для тестовых данных
    test_loaders = ChromosomeDataLoader[]
    @showprogress desc="Loading test genomes" for dir in test_directories
        d = digitize_genome_one_side(dir)
        all_dna = vcat(d.dna_pos, d.dna_neg)
        all_borders = vcat(d.borders_pos, d.borders_neg)
        for (dna, borders) in zip(all_dna, all_borders)
            loader = ChromosomeDataLoader(dna, borders, window_size, step)
            push!(test_loaders, loader)
        end
    end
    test_data = create_all_windows(test_loaders)
    test_dataloader = Flux.DataLoader(test_data, batchsize=batch_size, shuffle=false, collate=collate)

    # Создаем одну модель
    model = create_model()
    opt = ADAM(0.001)
    state = Flux.setup(opt, model)

    # Обучаем модель с оценкой на тестовых данных
    train_model!(model, state, train_dataloader, test_dataloader, epochs, weight_pos)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end