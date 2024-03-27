import pandas as pd
import numpy as np
from snakemake.io import expand
import pyranges as pr
import os
import cerberus

def assign_sector(df):
    df['sector'] = 'simple'

    df.loc[df.tss_ratio > 0.5, 'sector'] = 'tss'
    df.loc[df.tes_ratio > 0.5, 'sector'] = 'tes'
    df.loc[df.spl_ratio > 0.5, 'sector'] = 'splicing'

    # mixed genes
    df.loc[(df.sector=='simple')&(df.n_iso>1), 'sector'] = 'mixed'

    return df

def get_centroids(ca,
                  source='sample_det',
                  ver=None,
                  **kwargs):

    # subset on source
    df = ca.triplets.loc[ca.triplets.source == source].copy(deep=True)

    # limit only to columns that are important
    keep_cols = ['gname', 'gid',
                 'n_tss', 'n_tes', 'n_ic', 'splicing_ratio',
                 'tss_ratio', 'tes_ratio', 'spl_ratio',
                 'n_iso']
    df = df[keep_cols]

    # get centroid
    df = df.groupby(['gname', 'gid']).mean().reset_index()
    df = assign_sector(df)

    # add the centroids to the ca.triplet
    df['source'] = source+'_centroid'
    ca.triplets = pd.concat([ca.triplets, df], axis=0)

    return ca

def compute_dists(cas,
                  sources,
                  gene_merge=['gname', 'gid'],
                  rm_1_isos=[False, False],
                  gene_subsets=[None, None],
                  ver=[None, None]):
    """
    Compute the distance between source1 and source2.
    Also compute the Z-score.
    """

    def preproc_ca(ca,
                   source,
                   rm_1_iso,
                   gene_subset,
                   ver):
        """
        Preprocess cerberus annot according to input settings
        """

        # get triplets for source
        df = ca.triplets.loc[ca.triplets.source == source].copy(deep=True)

        # if requested, remove triplets w/ only 1 isoform
        if rm_1_iso:
            df = df.loc[df.n_iso > 1]

        # limit to target genes
        if gene_subset:
            gene_df, _, _ = get_gtf_info(how='gene',
                                         ver=ver,
                                         add_stable_gid=True)
            gene_df = gene_df[['gid_stable', 'biotype']]
            df = df.merge(gene_df, how='left',
                            left_on='gid', right_on='gid_stable')
            df = df.loc[df.biotype==gene_subset]

        return df

    if len(cas) > 2:
        print('Can only compute dists for 2 cerberus annots')
        return None

    dfs = []
    for i in range(len(cas)):
        dfs.append(preproc_ca(cas[i],
                               sources[i],
                               rm_1_isos[i],
                               gene_subsets[i],
                               ver[i]))





#     # get triplets for each source
#     df1 = ca.triplets.loc[ca.triplets.source == source1].copy(deep=True)
#     df2 = ca.triplets.loc[ca.triplets.source == source2].copy(deep=True)

#     # if requested, remove triplets w/ only one isoform
#     if rm_1_iso_1:
#         df1 = df1.loc[df1.n_iso > 1]
#     if rm_1_iso_2:
#         df2 = df2.loc[df2.n_iso > 1]

#     # limit to target genes
#     if gene_subset:
#         gene_df, _, _ = get_gtf_info(how='gene',
#                                      ver=ver,
#                                      add_stable_gid=True)
#         gene_df = gene_df[['gid_stable', 'biotype']]

#         # df1
#         df1 = df1.merge(gene_df, how='left',
#                         left_on='gid', right_on='gid_stable')
#         df1 = df1.loc[df1.biotype==gene_subset]

#         # df2
#         df2 = df2.merge(gene_df, how='left',
#                         left_on='gid', right_on='gid_stable')
#         df2 = df2.loc[df2.biotype==gene_subset]

    # merge dfs on gene info
    df = dfs[0].merge(dfs[1], how='inner',
                      on=gene_merge,
                      suffixes=(f'_{sources[0]}', f'_{sources[1]}'))

    # compute distances
    # pandarallel.initialize(nb_workers=8, verbose=1)
    df['dist'] = df.apply(simplex_dist,
                           args=(f'_{sources[0]}',
                                 f'_{sources[1]}'),
                           axis=1)
    df.dist = df.dist.fillna(0)

    # compute z_scores
    df['z_score'] = st.zscore(df.dist.tolist())

    return df



def get_dataset_df_col(wc, df, col):
    val = df.loc[df.dataset==wc.dataset, col].values[0]
    return val

def parse_input_config(config, config_file, subset):
    df = pd.read_csv(config_file, sep='\t')
    df['dataset_id'] = df.fname.str.split('.', n=1, expand=True)[0]
    df['subset'] = np.nan
    df.loc[df.fname.str.contains('all.gff.gz'), 'subset'] = 'all'
    df.loc[df.fname.str.contains('cagePolyASupported.gff.gz'), 'subset'] = 'supported'

    # add the metadata
    species = ['human', 'mouse']
    meta_df = pd.DataFrame()
    for s in species:
        temp = pd.read_csv(expand(config['data']['meta'], species=s)[0],
                           sep='\t')
        temp['species'] = s
        meta_df = pd.concat([meta_df, temp], axis=0)

    df = df.merge(meta_df,
                  left_on='dataset_id',
                  right_on='dataset',
                  how='left')
    df['dataset'] = df['tissue']+'_'+\
                 df['age']+'_'+\
                 df['platform']+'_'+\
                 df['capture']+'_'+\
                 df['species']+'_'+\
                 df['subset']

    df = df.loc[df.subset==subset]
    df['link'] = df.path+df.fname

    # assign a cerberus run to each dataset
    df['cerberus_run'] = df.sort_values(by=['species', 'dataset'],
                                 ascending=[True, True],
                                 na_position='last')\
                                 .groupby(['species'])\
                                 .cumcount() + 1

    df = df.sort_values(by='cerberus_run', ascending=True)
    return df

def get_ab_from_gff(gff_file, ofile):
    """
    Get abundance from the tagged values in the GFF file
    """
    df = pr.read_gff(gff_file).as_df()
    df = df[['transcript_id', 'flrpm', 'rpm']]
    df = df.drop_duplicates()
    assert len(df.loc[df.transcript_id.duplicated(keep=False)].index) == 0
    df.to_csv(ofile, sep='\t', index=False)

def format_tmerge_ab(ab, dataset, metric, ofile):
    """
    Format abundance file from tmerge GFF into a format
    that Cerberus can deal with
    """
    df = pd.read_csv(ab, sep='\t')
    df['annot_transcript_id'] = df['transcript_id']
    df['annot_transcript_name'] = df['transcript_id']
    df['transcript_ID'] = df['transcript_id']
    df.drop('transcript_id', axis=1, inplace=True)
    df[dataset] = df[metric]
    df.drop(['flrpm', 'rpm'], axis=1, inplace=True)
    df.to_csv(ofile, sep='\t', index=False)

def gff_rm_sirv(gff_file, ofile):
    """
    Remove SIRV chromosomes from GFF file
    """
    df = pr.read_gtf(gff_file).as_df()
    df = df.loc[df.Chromosome!='SIRVome_isoforms']
    df = pr.PyRanges(df)
    df.to_gtf(ofile)

def gff_fix_chr_names(gff_file, ofile, chr_map):
    """
    Fix the chromosome names of human gffs
    """
    m = pd.read_csv(chr_map, sep='\t')
    df = pr.read_gff(gff_file).as_df()
    df = df.merge(m, how='left',
                  left_on='Chromosome',
                  right_on='gff_chr')
    df.drop(['Chromosome', 'gff_chr'], axis=1, inplace=True)
    df.rename({'fa_chr': 'Chromosome'}, axis=1, inplace=True)
    df = pr.PyRanges(df)
    df.fillna(0, inplace=True)
    df.to_gtf(ofile)

def rm_multi_gene_ts(gff, ofile):
    """
    Remove transcripts from SQANTI GFF that matched with multiple genes.
    These are few in number (~10%)
    """
    df = pr.read_gff(gff).as_df()
    df['gene_count'] = df.gene_id.str.count('_')+1
    df = df.loc[df.gene_count == 1]
    df = pr.PyRanges(df)
    df.to_gtf(ofile)

def agg_cerb_abs(ab_files, ofile):
    """
    Put all abundance information from each sample in a
    transcripts x samples matrix
    """
    merge_cols = ['annot_transcript_name',
                  'annot_transcript_id']
    for i, f in enumerate(ab_files):
        temp = pd.read_csv(f, sep='\t')
        temp.drop('transcript_ID', axis=1, inplace=True)
        if i == 0:
            df = temp
        else:
            df = df.merge(temp, how='outer',
                          on=merge_cols)
    df.to_csv(ofile, sep='\t', index=False  )

def calculate_triplets(h5, ab, min_tpm, ofile, novel_genes=False):
    """
    Calculate triplets for genes based on expressed transcripts

    Parameters:
        h5 (str): Path to cerberus h5
        ab (str): Path to abundance matrix
        min_tpm (float): Min tpm for a transcript to be expressed
        ofile (str): Path to output cerberus h5
        novel_genes (bool): Keep novel genes. Default = Falses
    """
    # read in abundance
    df = pd.read_csv(ab, sep='\t')
    df = df.fillna(0)

    # read in cerberus annot
    ca = cerberus.read(h5)

    # remove novel genes?
    if not novel_genes:
        df = df.loc[~(df.annot_transcript_id.str.contains('novelGene'))]

    # determine whether each transcript is detected
    df.set_index(['annot_transcript_id', 'annot_transcript_name'], inplace=True)
    df = df >= min_tpm

    # loop through thing
    for c in df.columns.tolist():
        temp = df.loc[df[c] == True].copy(deep=True)
        temp.reset_index(inplace=True)
        tids = temp.annot_transcript_id.tolist()
        triplets = ca.get_subset_triplets(tids, c)
        ca.add_triplets(triplets)

    ca.write(ofile)
