import cerberus
tss = 'data/human/tss_agg.bed'
tes = 'data/human/tes_agg.bed'
ic = 'data/human/ic_agg.tsv'
h5 = 'test.h5'
cerberus.write_reference(tss,tes,ic,h5)
