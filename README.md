```bash
snakemake \
  -s snakemake/human_snakefile.smk \
  -j 10 \
  --latency-wait 120 \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END, --time=72:00:00" -n
```

Links:
* [CLS GTF](https://crgcnag-my.sharepoint.com/personal/gkaur_crg_es/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fgkaur%5Fcrg%5Fes%2FDocuments%2FAttachments%2FHv3%5FmasterTable%2Bgencodev42%5Flt1o%5Ftagged%5Floci%2Egtf%2Egz&parent=%2Fpersonal%2Fgkaur%5Fcrg%5Fes%2FDocuments%2FAttachments&ga=1)
