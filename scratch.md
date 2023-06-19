# rm sirvs
```python
import pyranges as pr
gtf='/dfs8/pub/freese/mortazavi_lab/data/230526_gencode/SIDMWbEPP0101.splicing_status-all.endSupport-all.gff'
df = pr.read_gff(gtf).as_df()
df = df.loc[df.Chromosome!='SIRVome_isoforms']
df = pr.PyRanges(df)
df.to_gtf('/dfs8/pub/freese/mortazavi_lab/data/230526_gencode/SIDMWbEPP0101_no_sirv.gtf')
```

<!-- # add transcript entries
```bash
# conda activate SQANTI3.env
# gffread \
#   -T /dfs8/pub/freese/mortazavi_lab/data/230526_gencode/SIDMWbEPP0101_no_sirv.gff \
#   -O /dfs8/pub/freese/mortazavi_lab/data/230526_gencode/SIDMWbEPP0101_no_sirv_t.gtf
conda activate SQANTI3.env
gffread \
  -E -T /dfs8/pub/freese/mortazavi_lab/data/230526_gencode/SIDMWbEPP0101_no_sirv.gtf \
  -o- > /dfs8/pub/freese/mortazavi_lab/data/230526_gencode/SIDMWbEPP0101_no_sirv_t.gtf
```
# try exogenously running gtf2genepred
```bash
/dfs8/pub/freese/mortazavi_lab/bin/SQANTI3/utilities/gtfToGenePred \
  /data/homezvol1/freese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.annotation.gtf \
  /dfs8/pub/freese/mortazavi_lab/data/230526_gencode/refAnnotation_/data/homezvol1/freese/mortazavi_lab/data/230526_gencode/test_sqanti.genePred', '-genePredExt', '-allErrors', '-ignoreGroupsWithoutExons'
``` -->

```bash
annot=/dfs8/pub/freese/mortazavi_lab/ref/gencode.vM21/gencode.vM21.annotation.gtf
ref=/dfs8/pub/freese/mortazavi_lab/ref/mm10/mm10.fa
# gtf=/dfs8/pub/freese/mortazavi_lab/data/230526_gencode/SIDMWbEPP0101_no_sirv_t.gtf
gtf=/dfs8/pub/freese/mortazavi_lab/data/230526_gencode/SIDMWbEPP0101_no_sirv.gtf
opref=test_sqanti
sqpath=/dfs8/pub/freese//mortazavi_lab/bin/SQANTI3/
conda activate SQANTI3.env
python ${sqpath}sqanti3_qc.py \
  $gtf \
   $annot \
   $ref \
   --force_id_ignore \
   --aligner_choice minimap2 \
   --skipORF \
   -o $opref
```

```bash
wget -O data/human/Placenta_Placenta_pacBioSII_post-capture_human_all.gff.gz https://public-docs.crg.es/rguigo/Data/gkaur/LyRic_CLS3_TM_perTissue/SIDHPlPPC0101.splicing_status-all.endSupport-all.gff.gz --user=user_cls --password=Gencode@CLS_2022
```
