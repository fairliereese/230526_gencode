rule cerberus_agg_ics_cfg:
    resources:
        threads = 1,
        nodes = 1
    run:
        refs = [params.ref for i in range(2)]
        df = pd.DataFrame()
        df['fname'] = [input.ref_ics, input.ics]
        df['ref'] = refs
        df['sources'] = params.sources
        df.to_csv(output.cfg, sep=',', header=None, index=False)

rule cerberus_agg_ics_cli:
    resources:
        threads = 2,
        nodes = 1
    shell:
        """
        cerberus agg_ics \
            --input {input.cfg} \
            -o {output.ics}
        """

use rule cerberus_agg_ics_cfg as cerb_agg_ics_cfg_lr with:
    input:
        ref_ics = lambda wc: get_prev_cerb_entry(wc, p_df,
                                                  config['analysis']['cerberus']['agg']['ics'],
                                                  config,
                                                  p_dir),
        ics = p_dir+config['cerberus']['ics']
    params:
        ref = False,
        sources = lambda wc:['cerberus', get_df_col(wc, df, 'source')]
    output:
        cfg = config['analysis']['cerberus']['agg']['ics_cfg']

use rule cerberus_agg_ics_cli as cerb_agg_ics_cfg_cli_lr with:
    input:
        ref_ics = lambda wc: get_prev_cerb_entry(wc, p_df,
                                                  config['analysis']['cerberus']['agg']['ics'],
                                                  config,
                                                  p_dir),
        ics = p_dir+config['cerberus']['ics'],
        cfg = config['analysis']['cerberus']['agg']['ics_cfg']
    output:
        ics = config['analysis']['cerberus']['agg']['ics']
