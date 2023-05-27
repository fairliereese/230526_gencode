import pandas as pd
import cerberus
import warnings
warnings.filterwarnings('ignore')

configfile: 'config.yml'

end_types = ['tss', 'tes']
config_fname = 'human_config.tsv'

df = pd.read_csv(config_fname, sep='\t')
species = df.species.tolist()[0]
source = df.source.tolist()[0] # would have to change in the future

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
    input:
           # expand(config['data']['cerberus']['ref'],
           #        zip,
           #        species=species)
           expand(config['data']['cerberus']['annot'],
                  zip,
                  species=species,
                  source=source)

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
        tss = lambda wc: expand(config['data']['cerberus']['agg_ends'],
                                zip,
                                species=wc.species,
                                end_type='tss')[0],
        tes = lambda wc: expand(config['data']['cerberus']['agg_ends'],
                                zip,
                                species=wc.species,
                                end_type='tes')[0],
        ic = config['data']['cerberus']['agg_ics']
    resources:
        mem_gb = 16,
        threads = 1
    output:
        h5 = config['data']['cerberus']['ref']
    run:
        cerberus.write_reference(input.tss,
                                 input.tes,
                                 input.ic,
                                 output.h5)

rule annot_ca:
    input:
        gtf = config['data']['gtf'],
        ref_h5 = config['data']['cerberus']['ref']
    # params:
    #     source = lambda wc:get_df_col(wc, df, 'source'),
    #     gene_source = lambda wc:get_df_col(wc, df, 'source'), # would have to change this mebbe
    resources:
        mem_gb = 32,
        threads = 4
    output:
        annot_h5 = config['data']['cerberus']['annot']
    run:
        cerberus.annotate_transcriptome(input.gtf,
                                        input.ref_h5,
                                        wildcards.source,
                                        None,
                                        output.annot_h5)
