# 0. Env
```shell
$ conda env create -n vkr -c bioconda treecluster ncbi-datasets-cli
$ conda activate vkr
```

# 1. Clustering
```shell
$ ./scripts/cluster_genomes.jl -i DATA/GTDB -o DATA/clusters
```


## Top 100 annotated:
```shell
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
DATA/GTDB/bac120_metadata_r220.tsv | \
sort -nr | \
head -n 100 > DATA/GTDB/top100_repeated_species.tsv
```

# 2. Downloading genomes
```shell
$ julia -t 10 scripts/download_genomes.jl -i DATA/cluseter/cg_hq.tsv -o DATA/genomes/ -e -p 10
```
