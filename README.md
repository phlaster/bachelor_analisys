# 1. Clustering
```bash
$ mkdir -p DATA/tree_GTDB

$ ./scripts/cluster_genomes.jl --InputDirectory DATA/tree_GTDB -O DATA/tree_GTDB/filtered_accessions -A
```
# 2. Downloading genomes
```bash
$ julia -t 10 scripts/download_genomes.jl -i DATA/tree_GTDB/filtered_accessions/filtered_accessions.txt -o DATA/genomes/ -e                                            
```

## Top 100 annotated:
```bash
$ awk -F"\t" '                                      
NR==1 {
  for(i=1;i<=NF;i++) {
    if($i=="gtdb_taxonomy") { tax_idx=i; break }
  }
  next
}
{
  if($tax_idx ~ /;s__[^;]+$/ || $tax_idx ~ /^s__[^;]+$/) {
    count[$tax_idx]++
  }
}
END {
  for(t in count)
    printf("%d\t%s\n", count[t], t)
}' \
DATA/tree_GTDB/bac120_metadata_r220.tsv | \
sort -nr | \
head -n 100 > top100_repeated_species.tsv
```

