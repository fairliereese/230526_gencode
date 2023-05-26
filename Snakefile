import pandas as pd
import cerberus
import warnings
warnings.filterwarnings('ignore')

configfile: 'snakemake/config.yml'

end_types = ['tss', 'tes']
config_fname = 'human_config.tsv'

df = pd.read_csv(config_fname, sep='\t')
species = df.species.tolist()[0]

if len(df.source.unique()) < len(df.source.tolist()):
    raise ValueError('Sources must have unique names')

def get_df_col(wc, df, col):
    if 'slack' in col:
        if wc.end_type == 'tss':
            col = 'tss_slack'
        elif wc.end_type == 'tes':
            col = 'tes_slack'
    elif 'extend' in col:
        if wc.end_type == 'tss':
            col = 'tss_extend'
        elif wc.end_type == 'tes':
            col = 'tes_extend'
    val = df.loc[df.source==wc.source, col].values[0]
    if 'extend' in col or 'slack' in col:
        val = int(val)
    return val

rule all:
    expand(config['data']['cerberus']['ref'],
           zip,
           species=species)

################################################################################
############################### Cerberus features ##############################
################################################################################

rule gtf_to_bed:
    input:
        gtf = config['data']['gtf']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        ends = config['data']['cerberus']['ends']
    params:
        slack = lambda wc:get_df_col(wc, df, 'slack'),
        dist = lambda wc:get_df_col(wc, df, 'extend')
    run:
        cerberus.gtf_to_bed(input.gtf,
                            wildcards.end_type,
                            output.ends,
                            dist=params.dist,
                            slack=params.slack)

rule gtf_to_ics:
    input:
        gtf = config['data']['gtf']
    resources:
        mem_gb = 64,
        threads = 1
    output:
        ics = config['data']['cerberus']['ics']
    run:
        cerberus.gtf_to_ics(input.gtf,
                            output.ics)

### TODO - make dummy agg_ends rule that just cps
rule agg_ends:
    input:
        ends = lambda wc: expand(config['data']['cerberus']['ends'],
                                 species=wc.species,
                                 source='cls',
                                 end_type=wc.end_type)
    resources:
        mem_gb = 8,
        threads = 1,
    output:
        ends = config['data']['cerberus']['agg_ends']
    shell:
        """
        cp {input.ends} {output.ends}
        """

rule agg_ics:
    input:
        ics = lambda wc: expand(config['data']['cerberus']['ics'],
                                 species=wc.species,
                                 source='cls')
    resources:
        mem_gb = 8,
        threads = 1,
    output:
        ics = config['data']['cerberus']['agg_ics']
    shell:
        """
        cp {input.ics} {output.ics}
        """

################################################################################
################################# Cerberus obj #################################
################################################################################

rule write_ca_ref:
    input:
        tss = lambda wc: expand(zip,
                                config['data']['cerberus']['agg_ends'],
                                species=wc.species,
                                end_type='tss'),
        tes = lambda wc: expand(zip,
                                config['data']['cerberus']['agg_ends'],
                                species=wc.species,
                                end_type='tes'),
        ic = config['data']['cerberus']['agg_ics']
    resources:
        mem_gb = 16,
        threads = 1
    output:
        h5 = config['data']['cerberus']['ref']
    run:
        write_reference(input.tss,
                        input.tes,
                        input.ic,
                        output.h5)
