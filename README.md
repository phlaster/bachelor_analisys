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

Calculating chromosomes:
```
find DATA/genomes/genomes -mindepth 1 -maxdepth 1 -type d | xargs -I{} -P 10 bash -c '
    dir="{}"                     
    file="$dir/$(basename "$dir").fna"
    count=$(grep -c ">" "$file")
    touch "$dir/n_chroms=$count"
'
```

# X. Experiment results:
## CDS starts
train 1000 pseudmnt
test 100 pseudmnt
### 1 `2025-04-08`
```
epochs=12
pad=30
window=61
gamma=3
lr=[0.005 x 3, 0.002 x 3, 0.0008 x 3, 0.00032 x 3]
loss=[0.00037685703, 0.00022302769, 0.00020763752, 0.00019436877, 0.00018928376, 0.00018650031, 0.00018145659, 0.00017990194, 0.00017909888, 0.00017715435, 0.00017665839, 0.00017631172]
CM:
┌────────┬───────────┬────────┬───────────┐
│ tr\pred│   Cl 1    │  Cl 2  │     Total │
├────────┼───────────┼────────┼───────────┤
│ Class1 │ 509146819 │  56100 │ 509202919 │
│ Class2 │    110345 │ 120642 │    230987 │
│  Total │ 509257164 │ 176742 │ 509433906 │
└────────┴───────────┴────────┴───────────┘
Class 1:
Prescision = 0.683
Recall     = 0.522
F1         = 0.592
FDR        = 0.317
Specificity= 1.0
Support    = 230987
```
### 2 `2025-04-09T22:15:58.046`
```
time=14686 seconds
epochs=10
pad=25
window=51
gamma=3
lr=[0.005, 0.003, 0.0018, 0.00108, 0.000648, 0.0003888, 0.00023328, 0.000139968, 8.39808e-5, 5.038848e-5]
loss=[0.001215547, 0.000233086, 0.00021530787, 0.0002071661, 0.0002007706, 0.00019636184, 0.00019344417, 0.00019183647, 0.0001908592, 0.00019024262]
CM:
┌────────┬───────────┬────────┐
│ tr\pred│   Cl 1    │  Cl 2  │
├────────┼───────────┼────────┤
│ Class1 │ 509156069 │  46850 │
│ Class2 │    124657 │ 106330 │
└────────┴───────────┴────────┘
Class 1:
Prescision = 0.694151
Recall     = 0.460329
F1         = 0.553561
FDR        = 0.305849
Specificity= 0.999908
Support    = 230987
```
### 3 `2025-04-09T22:20:01.208`
```
time=15231 seconds
epochs=10
pad=35
window=71
gamma=3
lr=[0.005, 0.003, 0.0018, 0.00108, 0.000648, 0.0003888, 0.00023328, 0.000139968, 8.39808e-5, 5.038848e-5]
loss=[0.00049170887, 0.00021514256, 0.00019736601, 0.00018857144, 0.00018282037, 0.00017939335, 0.00017721848, 0.00017586554, 0.00017515, 0.00017457388]
CM:
┌────────┬───────────┬────────┐
│ tr\pred│   Cl 1    │  Cl 2  │
├────────┼───────────┼────────┤
│ Class1 │ 509164058 │  38861 │
│ Class2 │    119174 │ 111813 │
└────────┴───────────┴────────┘
Class 1:
Prescision = 0.742086
Recall     = 0.484066
F1         = 0.585928
FDR        = 0.257914
Specificity= 0.999924
Support    = 230987
```