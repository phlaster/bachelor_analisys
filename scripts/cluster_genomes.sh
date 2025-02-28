#!/bin/bash

# Парсинг аргументов
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -d, --directory DIR       Output directory for tree and metadata files (required)"
    echo "  -t, --threshold THRESH    Clustering threshold for TreeCluster.py (required)"
    echo "  -c, --completeness CMP    Minimum completeness percentage (default: 95)"
    echo "  -C, --contamination CNT   Maximum contamination percentage (default: 5)"
    echo "  -m, --method METHOD       Clustering method [avg_clade|max_clade|...] (default: avg_clade)"
    echo "  -o, --output FILE         Output file name for selected accessions (default: DIR/selected_ncbi_accessions.txt)"
    echo "  -h, --help                Show this help message"
    exit 1
}

# Инициализация переменных
DIR=""
COMPLETENESS=95
CONTAMINATION=5
METHOD="avg_clade"
THRESHOLD=""
OUTPUT=""

# Парсинг аргументов
while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--directory)
            DIR="$2"
            shift 2
            ;;
        -t|--threshold)
            THRESHOLD="$2"
            shift 2
            ;;
        -c|--completeness)
            COMPLETENESS="$2"
            shift 2
            ;;
        -C|--contamination)
            CONTAMINATION="$2"
            shift 2
            ;;
        -m|--method)
            METHOD="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
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

# Проверка обязательных аргументов
if [ -z "$DIR" ] || [ -z "$THRESHOLD" ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Проверка доступности утилит
for cmd in awk wget gunzip sort join TreeCluster.py; do
    if ! command -v $cmd &> /dev/null; then
        echo "Error: $cmd not found. Please install it and try again."
        exit 1
    fi
done

# Создание временных файлов
TMP_HQ=$(mktemp)
TMP_SHQ=$(mktemp)
TMP_CLUSTER=$(mktemp)
TMP_SORT_CLUSTER=$(mktemp)
TMP_JOINED=$(mktemp)
TMP_SELECTED=$(mktemp)

# Удаление временных файлов при выходе
cleanup() {
    rm -f "$TMP_HQ" "$TMP_SHQ" "$TMP_CLUSTER" 
    rm -f "$TMP_SORT_CLUSTER" "$TMP_JOINED" "$TMP_SELECTED"
}
trap cleanup EXIT

# Создание основной директории
mkdir -p "$DIR"

# Скачивание файлов при необходимости
download_file() {
    local url=$1
    local dest=$2
    if [ ! -f "$dest" ]; then
        echo "Downloading ${url}..."
        wget -q -O - "$url" | gunzip > "$dest.tmp" && mv "$dest.tmp" "$dest"
    fi
}

METADATA_URL="https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_metadata_r220.tsv.gz"
TREE_URL="https://data.ace.uq.edu.au/public/gtdb/data/releases/release220/220.0/bac120_r220.tree.gz"

TREE_FILE="bac120_r220.nwk"
METADATA_FILE="bac120_metadata_r220.tsv"

download_file "$TREE_URL" "$DIR/$TREE_FILE"
download_file "$METADATA_URL" "$DIR/$METADATA_FILE"

# Фильтрация метаданных
echo "Filtering genomes (completeness >${COMPLETENESS}%, contamination <${CONTAMINATION}%)..."
awk -F'\t' -v c="$COMPLETENESS" -v C="$CONTAMINATION" \
    'NR>1 && $6 > c && $7 < C {print $1}' "$DIR/$METADATA_FILE" > "$TMP_HQ"
COUNT_HQ=$(wc -l < "$TMP_HQ")
echo "High-quality genomes: $COUNT_HQ"

# Кластеризация дерева
echo "Clustering tree with method=${METHOD}, threshold=${THRESHOLD}..."
TreeCluster.py -i "$DIR/$TREE_FILE" -o "$TMP_CLUSTER" -m "$METHOD" -t "$THRESHOLD"

# Обработка кластеров
tail -n +2 "$TMP_CLUSTER" | awk -F'\t' '{print $1 "\t" $2}' | sort -k1,1 > "$TMP_SORT_CLUSTER"

# Объединение с отфильтрованными геномами
sort "$TMP_HQ" > "$TMP_SHQ"
join -1 1 -2 1 "$TMP_SHQ" "$TMP_SORT_CLUSTER" > "$TMP_JOINED"
COUNT_JOINED=$(wc -l < "$TMP_JOINED")
echo "Genomes in clusters: $COUNT_JOINED"

# Выбор представителей
sort -k2,2n "$TMP_JOINED" | awk '!seen[$2]++ {print $1}' > "$TMP_SELECTED"
COUNT_SELECTED=$(wc -l < "$TMP_SELECTED")
echo "Selected representatives: $COUNT_SELECTED"

# Определение выходного файла
if [ -z "$OUTPUT" ]; then
    OUTPUT="$DIR/selected_ncbi_accessions.txt"
fi

# Создание директории для выходного файла
output_dir=$(dirname "$OUTPUT")
mkdir -p "$output_dir"

if [ -f "$OUTPUT" ]; then
    echo "Error: Output file $OUTPUT already exists. Use -o to specify another name."
    exit 1
fi

# Извлечение accession-ов
sed 's/.*\(GC[AF]_[0-9]\+\.[0-9]\+\).*/\1/' "$TMP_SELECTED" > "$OUTPUT"
echo "Results saved to: $OUTPUT"
echo "Total accessions: $(wc -l < "$OUTPUT")"