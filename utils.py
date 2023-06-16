import pandas as pd
import numpy as np
from snakemake.io import expand


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
    print(len(df.dataset.unique()))
    # print(df.head())
    return df
