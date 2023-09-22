import pandas as pd
import numpy as np
from snakemake.io import expand
import pyranges as pr
import os


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


def gff_rm_sirv(gff_file, ofile):
    """
    Remove SIRV chromosomes from GFF file
    """
    df = pr.read_gff(gff_file).as_df()
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
    df.to_gtf(ofile)
