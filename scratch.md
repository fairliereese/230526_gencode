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
```-->

# try exogenously running gtf2genepred
```bash
/dfs8/pub/freese/mortazavi_lab/bin/SQANTI3/utilities/gtfToGenePred \
  /share/crsp/lab/seyedam/freese/mortazavi_lab/data/230526_gencode/ref/mouse/annot.gtf \
  /share/crsp/lab/seyedam/freese/mortazavi_lab/data/230526_gencode/refAnnotation_data/mouse/WBlood_Adult_ont_pre-capture_mouse_all_sqanti.genePred \
  -genePredExt \
  -allErrors \
  -ignoreGroupsWithoutExons
```

```bash
python /dfs8/pub/freese//mortazavi_lab/bin/SQANTI3/sqanti3_qc.py \
  data/mouse/ESC_Embryo_pacBioSII_pre-capture_mouse_all_no_sirv.gff \
  ref/mouse/annot.gtf \
  ref/mouse/ref.fa \
  -d data/mouse/ESC_Embryo_pacBioSII_pre-capture_mouse_all_sqanti/ \
  --force_id_ignore \
  --aligner_choice minimap2 \
  --skipORF \
  -o data/mouse/ESC_Embryo_pacBioSII_pre-capture_mouse_all_sqanti

  python /dfs8/pub/freese//mortazavi_lab/bin/SQANTI3/sqanti3_qc.py \
    data/mouse/ESC_Embryo_pacBioSII_pre-capture_mouse_all_no_sirv.gff \
    ref/mouse/annot.gtf \
    ref/mouse/ref.fa \
    -d data/mouse/ESC_Embryo_pacBioSII_pre-capture_mouse_all/ \
    --force_id_ignore \
    --aligner_choice minimap2 \
    --skipORF \
    -o data/mouse/ESC_Embryo_pacBioSII_pre-capture_mouse_all_sqanti
```

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
