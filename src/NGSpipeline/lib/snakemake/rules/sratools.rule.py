rule fastq_dump:
    input:      
    output:     touch(config['wrk_folder'] + '/log/{experiment}.{run}.fastq_dump.done')
    params:     run_id = '{run}',
                output_folder = lambda ws: '{}/data/{}/'.format(config['wrk_folder'], ws.experiment)
    threads:    config['tools']['sra-tools']['threads']
    benchmark:  'log/{experiment}.{run}.fastq_dump.benchmark' 
    log:        outfile = 'log/{experiment}.{run}.fastq_dump.out',
                errfile = 'log/{experiment}.{run}.fastq_dump.err'
    conda:      config['conda_env_config']
    shell:      '''
                {config[tools][sra-tools][fastq-dump][call]}    \
                    -O {params.output_folder}                   \
                    {params.run_id}                             \
                    1> {log.outfile} 2> {log.errfile}
                '''
