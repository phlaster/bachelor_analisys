# 1. Clustering
```bash
$ mkdir -p DATA/tree_GTDB

$ ./scripts/cluster_genomes.sh -d DATA/tree_GTDB -c 99 -C 1 -t 0.5 -o DATA/tree_GTDB/selected_ncbi_accessions.txt
```
# 2. Downloading genomes
```bash
$ ./scripts/download_genomes.sh -i DATA/tree_GTDB/selected_ncbi_accessions.txt -o DATA/ -e
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

# 3. Quality check
```bash
$ mkdir -p DATA/checkm_input

$ find DATA/genomes/ncbi_dataset/data/ -name "*.fna" -exec ln -s {} DATA/checkm_input/ \;

$ checkm lineage_wf -t 32 -x fna DATA/checkm_input/ DATA/checkm_output/

$ checkm qa DATA/checkm_output/lineage.ms DATA/checkm_output -o 2 > DATA/checkm_output/quality_report.txt

$ checkm qa DATA/checkm_output/lineage.ms DATA/checkm_output -o 4 > DATA/checkm_output/detailed_report.txt
```

# 4. Tree
```bash
$ guppy tog DATA/checkm_output/storage/tree/concatenated.pplacer.json -o DATA/checkm_output/genome_tree_with_bins.nwk      
```





```
Genome assembly
    ASM430655v1
Taxon
    Rhizobium leguminosarum 
Strain
    SM52
Directory
    DATA/SM52
```