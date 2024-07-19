rule cerb_agg_ics_cfg:
  input:
      files = lambda wc:get_agg_settings(wc, 'file')
  resources:
      threads = 1,
      nodes = 1
  params:
      refs = False,
      sources = lambda wc:get_agg_settings(wc, 'source')
  output:
      cfg = config['data']['cerb']['agg_ics_cfg']
  run:
      df = pd.DataFrame()
      df['fname'] = input.files
      refs = [params.refs for i in range(len(input.files))]
      refs[0] = True
      df['ref'] = refs
      df['sources'] = params.sources
      df.to_csv(output.cfg, sep=',', header=None, index=False)

rule cerb_agg_ics:
    input:
        files = lambda wc:get_agg_settings(wc, 'file'),
        cfg = config['data']['cerb']['agg_ics_cfg']
    resources:
        threads = 16,
        nodes = 2
    output:
        ics = config['data']['cerb']['agg_ics']
    shell:
        """
        cerberus agg_ics \
            --input {input.cfg} \
            -o {output.ics}
        """
