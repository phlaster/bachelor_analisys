#!/usr/bin/env julia

const PROJECT_DIR = dirname(@__DIR__)
const UTILS_FILE = joinpath(PROJECT_DIR, "src", "utils.jl")

using Pkg
Pkg.activate(PROJECT_DIR, io=devnull)
include(UTILS_FILE)

using Flux
using LuxCUDA
using MLUtils
using Statistics

struct GenomeDataLoader
    dna_matrices::Vector{Matrix{UInt8}}  # Список матриц ДНК для каждой хромосомы
    starts::Vector{Vector{UInt8}}        # Список меток для каждой хромосомы
    window_size::Int                     # Размер окна
    step::Int                            # Шаг между окнами
end

function GenomeDataLoader(dna_matrix_pos::Vector{Matrix{UInt8}}, dna_matrix_neg::Vector{Matrix{UInt8}}, 
                         starts_pos::Vector{Vector{UInt8}}, starts_neg::Vector{Vector{UInt8}}, 
                         window_size::Int, step::Int)
    # Объединяем данные для прямой и обратной цепей, но сохраняем их как отдельные последовательности
    dna_matrices = vcat(dna_matrix_pos, dna_matrix_neg)
    starts = vcat(starts_pos, starts_neg)
    return GenomeDataLoader(dna_matrices, starts, window_size, step)
end

# Функция для получения окна с учётом кольцевой структуры
function get_window(dna_matrix::Matrix{UInt8}, starts::Vector{UInt8}, i::Int, window_size::Int, seq_len::Int)
    if seq_len >= window_size
        indices = [mod1(j, seq_len) for j in i:(i + window_size - 1)]
        window_dna = dna_matrix[indices, :]
        window_starts = starts[indices]
    else
        repeats = ceil(Int, window_size / seq_len)
        extended_dna = repeat(dna_matrix, repeats, 1)
        extended_starts = repeat(starts, repeats)
        window_dna = extended_dna[i:(i + window_size - 1), :]
        window_starts = extended_starts[i:(i + window_size - 1)]
    end
    return window_dna, window_starts
end

# Преобразование данных в формат Flux
function to_flux_format(window_dna::Matrix{UInt8}, window_starts::Vector{UInt8})
    X = Float32.(window_dna)
    X = permutedims(X, (2, 1))
    y = Float32.(window_starts)
    return X, y
end

function collate(batch)
    X = cat([b[1] for b in batch]..., dims=3)
    X = permutedims(X, (2, 1, 3))
    y = cat([b[2] for b in batch]..., dims=2)
    return X, y
end

function create_dataloader(loader::GenomeDataLoader, batch_size::Int)
    data = Tuple{Matrix{Float32}, Vector{Float32}}[]
    for (dna_matrix, starts) in zip(loader.dna_matrices, loader.starts)
        seq_len = size(dna_matrix, 1)
        for pos in 1:loader.step:seq_len
            window_dna, window_starts = get_window(dna_matrix, starts, pos, loader.window_size, seq_len)
            X, y = to_flux_format(window_dna, window_starts)
            push!(data, (X, y))
        end
    end
    return Flux.DataLoader(data, batchsize=batch_size, shuffle=true, collate=collate)
end

function Base.iterate(loader::GenomeDataLoader, state=nothing)
    if state === nothing
        state = (1, 1)  # (seq_index, pos_index)
    end
    seq_index, pos_index = state
    
    # Проверяем, закончились ли последовательности
    if seq_index > length(loader.dna_matrices)
        return nothing
    end
    
    dna_matrix = loader.dna_matrices[seq_index]
    starts = loader.starts[seq_index]
    seq_len = size(dna_matrix, 1)
    
    # Получаем окно
    window_dna, window_starts = get_window(dna_matrix, starts, pos_index, loader.window_size, seq_len)
    X, y = to_flux_format(window_dna, window_starts)
    
    # Следующая позиция
    next_pos = pos_index + loader.step
    if next_pos + loader.window_size - 1 > seq_len
        next_seq = seq_index + 1
        next_pos = 1
    else
        next_seq = seq_index
    end
    next_state = (next_seq, next_pos)
    
    return (X, y), next_state
end

Base.length(loader::GenomeDataLoader) = sum(seq -> cld(size(seq, 1), loader.step), loader.dna_matrices)
Base.eltype(loader::GenomeDataLoader) = Tuple{Matrix{Float32}, Vector{Float32}}


function check_data_consistency(loader::GenomeDataLoader)
    for (dna, starts) in zip(loader.dna_matrices, loader.starts)
        @assert size(dna, 1) == length(starts) "Mismatch between DNA matrix rows and starts vector length"
        @assert size(dna, 2) == 4 "DNA matrix must have 4 columns for one-hot encoding"
    end
    println("Data consistency check passed")
end



function train_model!(model, dataloader, epochs=10)
    for epoch in 1:epochs
        for (X, y) in dataloader
            X_gpu = X |> gpu
            y_gpu = y |> gpu
            gs = gradient(model -> loss(model(X_gpu), y_gpu), model)
            Flux.update!(state, model, gs[1])
        end
        println("Эпоха $epoch завершена")
    end
end


genome_directories = "DATA/genomes/genomes/" .* readdir("DATA/genomes/genomes")[1:10]


all_dna_matrix_pos = Vector{Matrix{UInt8}}[]
all_dna_matrix_neg = Vector{Matrix{UInt8}}[]
all_starts_pos = Vector{Vector{UInt8}}[]
all_ends_pos = Vector{Vector{UInt8}}[]
all_starts_neg = Vector{Vector{UInt8}}[]
all_ends_neg = Vector{Vector{UInt8}}[]


for dir in genome_directories
    dna_matrix_pos, dna_matrix_neg, starts_pos, ends_pos, starts_neg, ends_neg = digitize_genome(dir)
    
    push!(all_dna_matrix_pos, dna_matrix_pos)
    push!(all_dna_matrix_neg, dna_matrix_neg)
    push!(all_starts_pos, starts_pos)
    push!(all_ends_pos, ends_pos)
    push!(all_starts_neg, starts_neg)
    push!(all_ends_neg, ends_neg)
end

all_dna_matrices = vcat(all_dna_matrix_pos..., all_dna_matrix_neg...)
all_starts = vcat(all_starts_pos..., all_starts_neg...)

# Параметры загрузчика
window_size = 200
step = 50
batch_size = 32

loader = GenomeDataLoader(all_dna_matrices, all_starts, window_size, step)

check_data_consistency(loader)

dataloader = create_dataloader(loader, batch_size)

model = Chain(
    Conv((5,), 4 => 32, relu, pad=SamePad()),  # Ожидает (window_size, 4, batch_size)
    Conv((5,), 32 => 64, relu, pad=SamePad()),
    Conv((1,), 64 => 1, identity, pad=SamePad()),
    x -> sigmoid.(x),
    x -> dropdims(x, dims=2)  # Удаляем размерность канала: (window_size, batch_size)
) |> gpu

state = Flux.setup(opt, model)

loss(y_pred, y) = mean(Flux.binarycrossentropy.(y_pred, y))

opt = ADAM(0.001)  # Learning rate = 0.001

train_model!(model, dataloader, 10)