################################################################################
####################### Getting transcript / gene info #########################
################################################################################

def get_polya_cats():
    return ['protein_coding', 'lncRNA', 'pseudogene']

def get_biotype_map():
    """
    Get a dictionary mapping each gene type to a more general biotype
    """
    map = {'protein_coding': ['protein_coding'],
           'lncRNA': ['lincRNA',
                      'lncRNA',
                      'processed_transcript',
                      'sense_intronic',
                      '3prime_overlapping_ncRNA',
                      'bidirectional_promoter_lncRNA',
                      'sense_overlapping',
                      'non_coding',
                      'macro_lncRNA',
                      'antisense'],
           'pseudogene': ['unprocessed_pseudogene',
                          'translated_unprocessed_pseudogene',
                          'transcribed_unprocessed_pseudogene',
                          'processed_pseudogene',
                          'transcribed_processed_pseudogene',
                          'transcribed_unitary_pseudogene',
                          'unitary_pseudogene',
                          'polymorphic_pseudogene',
                          'pseudogene',
                          'translated_processed_pseudogene'],
           'miRNA': ['miRNA'],
           'other': ['snRNA', 'vault_RNA',
                     'misc_RNA', 'TEC',
                     'snoRNA', 'scaRNA',
                     'rRNA_pseudogene', 'rRNA',
                     'IG_V_pseudogene',
                     'scRNA', 'IG_V_gene',
                     'IG_C_gene', 'IG_J_gene',
                     'sRNA', 'ribozyme',
                     'vaultRNA', 'TR_C_gene',
                     'TR_J_gene', 'TR_V_gene',
                     'TR_V_pseudogene', 'TR_D_gene',
                     'IG_C_pseudogene', 'IG_D_gene',
                     'IG_pseudogene', 'Mt_tRNA',
                     'Mt_rRNA', 'TR_J_pseudogene',
                     'IG_J_pseudogene']}
    return map

def get_transcript_info(gtf, tf_file, o):
    """
    Get a file with relevant information about each transcript in a gtf

    Parameters:
        gtf (str): Path to gtf
        o (str): Output file name
    """

    df = pr.read_gtf(gtf, as_df=True, duplicate_attr=True)

    # remove sirvs and erccs
    print(len(df.index))
    df = df.loc[(~df.Chromosome.str.contains('SIRV'))&~(df.Chromosome.str.contains('ERCC'))]
    print(len(df.index))

    # mane status
    mane_df = df.loc[df.Feature == 'transcript'].copy(deep=True)
    mane_df.tag.fillna('', inplace=True)
    mane_df['MANE_Select'] = mane_df.tag.str.contains('MANE_Select')
    mane_df['MANE_Plus_Clinical'] = mane_df.tag.str.contains('MANE_Plus_Clinical')
    mane_df = mane_df[['transcript_id', 'MANE_Select', 'MANE_Plus_Clinical']]
    mane_df.rename({'transcript_id':'tid'}, axis=1, inplace=True)

    # only exons
    df = df.loc[df.Feature == 'exon'].copy(deep=True)

    # readthrough transcripts
    df['readthrough_transcript'] = df.tag.str.contains('readthrough_transcript')

    # rename some columns
    m = {'gene_id': 'gid',
         'gene_name': 'gname',
         'transcript_id': 'tid',
         'gene_type': 'biotype'}
    df.rename(m, axis=1, inplace=True)

    map = get_biotype_map()

    beeps = []
    for key, item in map.items():
        beeps += item

    set(df.biotype.unique().tolist())-set(beeps)

    # pivot map
    biotype_map = {}
    for key, biotypes in map.items():
        for biotype in biotypes:
            biotype_map[biotype] = key

    # then add map to df
    df['biotype_category'] = df.biotype.map(biotype_map)

    df['exon_len'] = (df.Start-df.End).abs()+1

    cols = ['gid', 'gname', 'tid', 'exon_len', 'biotype', 'biotype_category', 'readthrough_transcript']
    df = df[cols]
    df_copy = df[['gid', 'gname', 'tid', 'biotype', 'biotype_category', 'readthrough_transcript']].copy(deep=True)
    df_copy = df_copy.drop_duplicates(keep='first')

    df = df.groupby('tid').sum().reset_index()
    df.rename({'exon_len': 't_len'}, axis=1, inplace=True)
    df = df.merge(df_copy, on='tid', how='left')

    # merge mane info
    df = df.merge(mane_df, how='left', on='tid')

    # # add TF info
    # df['tf'] = False
    # tf_df = pd.read_csv(tf_file, sep='\t')
    # tf_gids = tf_df['Gene stable ID'].unique().tolist()
    # df['gid_stable'] = df['gid'].str.split('.', expand=True)[0]
    # df.loc[df.gid_stable.isin(tf_gids), 'tf'] = True
    # df.drop('gid_stable', axis=1, inplace=True)

    # and save
    df.to_csv(o, sep='\t', index=False)

def get_gene_info(gtf, o):
    df = pr.read_gtf(gtf, as_df=True, duplicate_attr=True)

    # remove sirvs and erccs
    print(len(df.index))
    df = df.loc[(~df.Chromosome.str.contains('SIRV'))&~(df.Chromosome.str.contains('ERCC'))]
    print(len(df.index))

    # mane / mane clinical status
    df['MANE_Select'] = df.tag.str.contains('MANE_Select')
    df['MANE_Plus_Clinical'] = df.tag.str.contains('MANE_Plus_Clinical')
    temp = df[['gene_id', 'MANE_Select', 'MANE_Plus_Clinical']].copy(deep=True)
    df.drop(['MANE_Select', 'MANE_Plus_Clinical'], axis=1, inplace=True)
    temp = temp.groupby('gene_id').max().reset_index()

    # only gene, merge in that stuff
    df = df.loc[df.Feature == 'gene'].copy(deep=True)
    df = df.merge(temp, how='left', on='gene_id')

    # rename some columns
    m = {'gene_id': 'gid',
         'gene_name': 'gname',
         'transcript_id': 'tid',
         'gene_type': 'biotype'}
    df.rename(m, axis=1, inplace=True)

    map = get_biotype_map()

    beeps = []
    for key, item in map.items():
        beeps += item

    set(df.biotype.unique().tolist())-set(beeps)

    # pivot map
    biotype_map = {}
    for key, biotypes in map.items():
        for biotype in biotypes:
            biotype_map[biotype] = key

    # then add map to df
    df['biotype_category'] = df.biotype.map(biotype_map)

    # gene length
    df['length'] = (df.Start-df.End).abs()

    # # add TF info
    # df['tf'] = False
    # tf_df = pd.read_csv(tf_file, sep='\t')
    # tf_gids = tf_df['Gene stable ID'].unique().tolist()
    # df['gid_stable'] = df['gid'].str.split('.', expand=True)[0]
    # df.loc[df.gid_stable.isin(tf_gids), 'tf'] = True
    # df.drop('gid_stable', axis=1, inplace=True)

    # and save
    df = df[['gid', 'gname', 'length', 'biotype', 'biotype_category', 'MANE_Select', 'MANE_Plus_Clinical']]
    df.to_csv(o, sep='\t', index=False)

rule get_t_info:
    resources:
        nodes = 1,
        threads = 1
    run:
        get_transcript_info(input.gtf,
                            input.tf_file,
                            output.o)

rule get_g_info:
    resources:
        nodes = 1,
        threads = 1
    run:
        get_gene_info(input.gtf,
                      output.o)

use rule get_g_info as gc_g_info_ref with:
  input:
      gtf = config['ref']['annot']
  output:
      o = config['ref']['g_info']

rule all_refs:
    input:
        expand(config['ref']['g_info'],
            species=df.species.tolist())
