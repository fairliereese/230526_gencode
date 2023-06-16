import pandas as pd
import cerberus
import warnings
from utils import *

warnings.filterwarnings('ignore')

configfile: 'config.yml'

end_types = ['tss', 'tes']

config_fname = '230614_config.tsv'
df = parse_input_config(config, config_fname, 'all')
datasets = df.dataset.tolist()
species = df.species.tolist()

rule all:
    input:
        expand(config['data']['gff_gz'],
               zip,
               species=species,
               dataset=datasets)

################################################################################
########################### Ref. processing ####################################
################################################################################

rule dl:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "wget -O {output.out} {params.link}"

rule dl_pass:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "wget -O {output.out} {params.link} --user={params.user} --password={params.pwd}"

rule gunzip:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "gunzip -c {input.gz} > {output.out}"

use rule dl_pass as dl_gff with:
    params:
        link = lambda wc: get_dataset_df_col(wc, df, 'link'),
        user = 'user_cls',
        pwd = 'Gencode@CLS_2022'
    output:
        out = config['data']['gff_gz']
