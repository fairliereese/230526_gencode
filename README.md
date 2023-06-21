```bash
snakemake \
  -s Snakefile \
  -j 100 \
  --latency-wait 120 \
  --use-conda \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n
```

```bash
conda activate base_clone
snakemake \
  -s Snakefile \
  -j 100 \
  --latency-wait 120 \
  -n
```

Links:
* [Human metadata](https://github.com/guigolab/gencode-cls-master-table/releases/latest/download/Hv3_metadata.tsv.gz)
* [Mouse metadata](https://github.com/guigolab/gencode-cls-master-table/releases/latest/download/Mv2_metadata.tsv.gz)

* For the metadata files, I manually added column names
