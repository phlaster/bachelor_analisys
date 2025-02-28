#!/bin/bash

# Парсинг аргументов
usage() {
    echo "Usage: $0 -i INPUT_FILE -o OUTPUT_DIR [-c CHUNK_SIZE] [-r MAX_RETRIES] [-e]"
    echo "Options:"
    echo "  -i, --input-file    Path to accessions file (required)"
    echo "  -o, --output-dir    Output directory (required)"
    echo "  -c, --chunk-size    Number of accessions per chunk (default: 50)"
    echo "  -r, --max-retries   Maximum download attempts per chunk (default: 3)"
    echo "  -e, --extract       Extract downloaded archives (default: false)"
    exit 1
}

# Инициализация переменных
INPUT_FILE=""
OUTPUT_DIR=""
CHUNK_SIZE=50
MAX_RETRIES=3
EXTRACT=0

# Проверка зависимостей
check_dependencies() {
    local deps=("split" "datasets" "unzip" "sed" "basename" "date")
    for cmd in "${deps[@]}"; do
        if ! command -v "$cmd" &> /dev/null; then
            echo "Error: Required command '$cmd' not found"
            exit 1
        fi
    done
}

# Парсинг аргументов командной строки
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input-file)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -c|--chunk-size)
            CHUNK_SIZE="$2"
            shift 2
            ;;
        -r|--max-retries)
            MAX_RETRIES="$2"
            shift 2
            ;;
        -e|--extract)
            EXTRACT=1
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Проверка обязательных параметров
if [ -z "$INPUT_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Input file and output directory are required"
    usage
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file not found: $INPUT_FILE"
    exit 1
fi

check_dependencies

# Создание директорий
mkdir -p "$OUTPUT_DIR/accessions_chunks"
mkdir -p "$OUTPUT_DIR/archives"
[ $EXTRACT -eq 1 ] && mkdir -p "$OUTPUT_DIR/genomes"

# Разбиение на чанки
echo "Splitting input file into chunks of $CHUNK_SIZE..."
split -l "$CHUNK_SIZE" -d "$INPUT_FILE" "$OUTPUT_DIR/accessions_chunks/accessions_chunk_"

# Переменная для отслеживания минимального времени загрузки
MIN_DURATION=-1

# Скачивание чанков
process_chunk() {
    local chunk_path="$1"
    local chunk_num=$(basename "$chunk_path" | sed 's/accessions_chunk_//')
    local archive_path="$OUTPUT_DIR/archives/genomes_chunk_${chunk_num}.zip"
    
    local attempts=0
    local success=0
    local duration=0

    while [ $attempts -lt $MAX_RETRIES ]; do
        # Удаление предыдущего архива
        rm -f "$archive_path"
        
        # Замер времени скачивания
        local start=$(date +%s)
        
        if datasets download genome accession \
            --inputfile "$chunk_path" \
            --include genome,gff3 \
            --filename "$archive_path" 2>/dev/null; then
            
            local end=$(date +%s)
            duration=$((end - start))
            
            # Проверка времени выполнения
            if [ $MIN_DURATION -ne -1 ] && [ $duration -gt $((MIN_DURATION * 10)) ]; then
                echo "Chunk $chunk_num: Timeout (${duration}s), retrying..."
                ((attempts++))
                continue
            fi
            
            # Обновление минимального времени
            if [ $MIN_DURATION -eq -1 ] || [ $duration -lt $MIN_DURATION ]; then
                MIN_DURATION=$duration
            fi
            
            success=1
            break
        else
            ((attempts++))
            echo "Chunk $chunk_num: Attempt $attempts failed"
        fi
    done

    if [ $success -eq 1 ]; then
        echo "Chunk $chunk_num: Downloaded successfully in ${duration}s"
    else
        echo "Chunk $chunk_num: Failed after $MAX_RETRIES attempts"
    fi
}

# Обработка всех чанков
for chunk in "$OUTPUT_DIR/accessions_chunks/accessions_chunk_"*; do
    process_chunk "$chunk"
done

# Распаковка архивов
if [ $EXTRACT -eq 1 ]; then
    echo "Extracting archives..."
    for zip in "$OUTPUT_DIR/archives/"*.zip; do
        unzip -o "$zip" -d "$OUTPUT_DIR/genomes" >/dev/null
    done
    echo "Extraction completed"
fi

echo "Download process finished successfully"