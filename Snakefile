import pandas as pd
import cerberus
import warnings
from utils import *

warnings.filterwarnings('ignore')

configfile: 'config.yml'

end_types = ['tss', 'tes']

# data config
config_fname = '230614_config.tsv'
df = parse_input_config(config, config_fname, 'all')

# todo
df = df.loc[df.species == 'mouse']

datasets = df.dataset.tolist()
species = df.species.tolist()
max_cerberus_run = len(datasets) # TODO modify to work w/ multiple species

cerb_tsv = 'cerberus.tsv'
cerb_settings = pd.read_csv(cerb_tsv, sep='\t')

def get_df_col(df, q_col, q_col_val, t_col):
    return df.loc[df[q_col]==q_col_val, t_col].tolist()

wildcard_constraints:
    dataset= '|'.join([re.escape(x) for x in datasets]),

rule all:
    input:
        # 'beep_env.out'
        # expand(expand(config['data']['cerb']['ends'],
        #     zip,
        #     species=species,
        #     dataset=datasets,
        #     allow_missing=True),
        #     end_mode=end_types),
        # expand(config['data']['cerb']['ics'],
        #         zip,
        #         species=species,
        #         dataset=datasets)
        # expand(config['data']['cerb']['ca_ref'],
        #        zip,
        #        species=species)
        expand(config['data']['ca_annot'],
               zip,
               species=species,
               dataset=dataset[-1],
               cerberus_run=max_cerberus_run)

# rule debug_envs:
#     conda:
#         "SQANTI3.env"
#     resources:
#         threads = 1,
#         mem_gb = 4
#     output:
#         out = 'beep_env.out'
#     shell:
#         """
#             python test.py > {output.out}
#         """

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

use rule gunzip as gunzip_gff with:
    input:
        gz = config['data']['gff_gz']
    output:
        out = config['data']['gff']

use rule dl as dl_annot with:
    params:
        link = lambda wc: config['ref'][wc.species]['annot_link']
    output:
        out = config['ref']['annot_gz']

use rule dl as dl_fa with:
    params:
        link = lambda wc: config['ref'][wc.species]['fa_link']
    output:
        out = config['ref']['fa_gz']

use rule gunzip as gunzip_annot with:
    input:
        gz = config['ref']['annot_gz']
    output:
        out = config['ref']['annot']

use rule gunzip as gunzip_fa with:
    input:
        gz = config['ref']['fa_gz']
    output:
        out = config['ref']['fa']

################################################################################
########################## Cerberus ref stuff ##################################
################################################################################

def get_cerb_settings(wc, df, col):
    """
    Get Cerberus end-related settings
    """
    col = f'{wc.end_mode}_{col}'
    return df[col].tolist()[0]

rule cerb_gtf_to_bed:
    resources:
        mem_gb = 64,
        threads = 1
    run:
        cerberus.gtf_to_bed(input.gtf,
                            wildcards.end_mode,
                            output.ends,
                            dist=params.dist,
                            slack=params.slack)

rule cerb_gtf_to_ics:
    resources:
        mem_gb = 64,
        threads = 1
    run:
        cerberus.gtf_to_ics(input.gtf,
                            output.ics)

use rule cerb_gtf_to_bed as ref_cerb_gtf_to_bed with:
    input:
        gtf = config['ref']['annot']
    output:
        ends = config['ref']['ends']
    params:
        slack = lambda wc:get_cerb_settings(wc, cerb_settings, 'slack'),
        dist = lambda wc:get_cerb_settings(wc, cerb_settings, 'dist')

use rule cerb_gtf_to_ics as ref_cerb_gtf_to_ics with:
    input:
        gtf = config['ref']['annot']
    output:
        ics = config['ref']['ics']

################################################################################
###################### Data formatting / cleaning ##############################
################################################################################

rule get_gff_ab:
    input:
        gff = config['data']['gff']
    resources:
        mem_gb = 64,
        threads = 4
    output:
        ab = config['data']['ab']
    run:
        get_ab_from_gff(input.gff, output.ab)

rule rm_sirv:
    input:
        gff = config['data']['gff']
    resources:
        mem_gb = 64,
        threads = 4
    output:
        gtf = config['data']['gtf_no_sirv']
    run:
        gff_rm_sirv(input.gff, output.gtf)

################################################################################
################################# SQANTI #######################################
################################################################################

rule sqanti:
    input:
        gtf = config['data']['gtf_no_sirv'],
        fa = config['ref']['fa'],
        annot = config['ref']['annot']
    resources:
        mem_gb = 128,
        threads = 8
    conda:
        "SQANTI3.env"
    params:
        sq_path = config['sqanti_path'],
        c_path = config['cupcake_path'],
        opref = config['data']['sqanti_gff'].split('/')[-1].split('_corrected')[0],
        odir = '/'.join(config['data']['sqanti_gff'].split('/')[:-1])+'/'
    output:
        gff = config['data']['sqanti_gff']
    shell:
        """
        mkdir -p {params.odir}
        python {params.sq_path}sqanti3_qc.py \
            {input.gtf} \
            {input.annot} \
            {input.fa} \
            -d {params.odir} \
            --report skip \
            --force_id_ignore \
            --aligner_choice minimap2 \
            --skipORF \
            --cupcake_path {params.c_path} \
            -o {params.opref}
        """

################################################################################
############################## Cerberus ########################################
################################################################################

use rule cerb_gtf_to_bed as data_cerb_gtf_to_bed with:
    input:
        gtf = config['data']['sqanti_gff']
    output:
        ends = config['data']['cerb']['ends']
    params:
        slack = lambda wc:get_cerb_settings(wc, cerb_settings, 'slack'),
        dist = lambda wc:get_cerb_settings(wc, cerb_settings, 'dist')

use rule cerb_gtf_to_ics as data_cerb_gtf_to_ics with:
    input:
        gtf = config['data']['sqanti_gff']
    output:
        ics = config['data']['cerb']['ics']

def get_agg_settings(wc, param='files'):
    settings = {}
    files = []
    sources = []

    if 'end_mode' in wc.keys():
        ic = False
    else:
        ic = True

    # get reference ics / tsss / tess first
    # first source is 'cerberus' to indicate
    # preservation of source / novelty
    if ic == True:
        files += [config['ref']['ics']]
    else:
        files += expand(config['ref']['ends'],
                        zip,
                        end_mode=wc.end_mode,
                        species=wc.species)
    sources += ['gencode']

    # then get the files from each study
    # use the datasets as the sources
    if ic == True:
        files += expand(expand(config['data']['cerb']['ics'],
               zip,
               dataset=get_df_col(df, 'species', wc.species, 'dataset'),
               allow_missing=True),
               species=wc.species)
    else:
        files += expand(expand(config['data']['cerb']['ends'],
               zip,
               dataset=get_df_col(df, 'species', wc.species, 'dataset'),
               allow_missing=True),
               species=wc.species,
               end_mode=wc.end_mode)
    sources += datasets

    settings['file'] = files
    settings['source'] = sources

    return settings[param]

rule cerb_agg_ends:
    input:
        files = lambda wc:get_agg_settings(wc, 'file')
    resources:
      threads = 4,
      mem_gb = 32
    params:
        add_ends = True,
        refs = False,
        slack = lambda wc:get_cerb_settings(wc, cerb_settings, 'agg_slack'),
        sources = lambda wc:get_agg_settings(wc, 'source')
    output:
        ends = config['data']['cerb']['agg_ends']
    run:
        refs = [params.refs for i in range(len(input.files))]
        add_ends = [params.add_ends for i in range(len(input.files))]
        cerberus.agg_ends(input.files,
                          add_ends,
                          refs,
                          params.sources,
                          wildcards.end_mode,
                          params.slack,
                          output.ends)

rule cerb_agg_ics:
  input:
      files = lambda wc:get_agg_settings(wc, 'file')
  resources:
    threads = 4,
    mem_gb = 32
  params:
      refs = False,
      sources = lambda wc:get_agg_settings(wc, 'source')
  output:
      ics = config['data']['cerb']['agg_ics']
  run:
      refs = [params.refs for i in range(len(input.files))]
      refs[0] = True
      cerberus.agg_ics(input.files,
                        refs,
                        params.sources,
                        output.ics)

rule cerb_write_ref:
    input:
        ic = config['data']['cerb']['agg_ics'],
        tss = lambda wc:expand(config['data']['cerb']['agg_ends'],
                               species=wc.species,
                               end_mode='tss')[0],
        tes = lambda wc:expand(config['data']['cerb']['agg_ends'],
                              species=wc.species,
                              end_mode='tes')[0]
    resources:
        threads = 4,
        mem_gb = 64
    output:
        h5 = config['data']['cerb']['ca_ref']
    run:
        cerberus.write_reference(input.tss,
                                 input.tes,
                                 input.ic,
                                 output.h5)

 ################################################################################
 ######################### Cerberus annot + ID replacement ######################
 ################################################################################

rule cerb_annot:
    resources:
        mem_gb = 64,
        threads = 16
run:
    cerberus.annotate_transcriptome(input.gtf,
                                 input.h5,
                                 params.source,
                                 params.gene_source,
                                 output.h5)

use rule cerb_annot as ref_cerb_annot with:
    input:
        h5 = config['data']['cerb']['ca_ref'],
        gtf = lambda wc: config['ref']['annot']
    params:
        source = 'gencode',
        gene_source = None
    output:
        h5 = config['ref']['ca_annot']

def get_prev_ca_annot(wc):
    # TODO - modify to work w/ multiple species
    if wc.cerberus_run == 0:
        ca = expand(config['ref']['ca_annot'],
                    zip,
                    species=wc.species)
    else:
        prev_run = int(wc.cerberus_run)-1
        prev_dataset = datasets[prev_run]
        ca = expand(config['data']['cerb']['ca_annot'],
                    zip,
                    dataset=prev_dataset,
                    cerberus_run=prev_run,
                    species=wc.species)
    return ca

# TODO need this to modify the same cerberus obj sequentially
use rule cerb_annot as study_cerb_annot with:
    input:
        h5 = lambda wc: get_prev_ca_annot(wc, cerberus_run)
        gtf = config['data']['gtf_no_sirv']
    params:
        source = lambda wc:wc.dataset,
        gene_source = 'gencode'
    output:
        h5 = config['data']['ca_annot']
