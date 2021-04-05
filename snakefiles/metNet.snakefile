rule metNet_run:
    input:
        input = config['']['input'],
        ouptut = config['']['output'],
        model = config['']['model'],
        classifications = config['']['classifications'],
        state = config['']['state'],

    threads: config['args']['max_threads_per_rule']

    conda: config['envs']['metNet']

    params:
        sample_name = config['args']['input'],
        config = config['config']['metNet']['files'][0],
        
        
    
    output:
        config['outdir']['metNet']+"output_directory/classifications/Out.txt",
        config['outdir']['metNet']+"output_directory/classifications/CM.png",
        config['outdir']['metNet']+"output_directory/classifications/PRHM.png",
    
    script:
        config['args']['mcc_path']+"/scripts/metNet.py"
