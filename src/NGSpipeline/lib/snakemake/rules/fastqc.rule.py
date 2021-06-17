def dump_or_gunzip(ws):
    if (ws.state == 'unfiltered'):
        if (config['experiments'][ws.experiment]['runs'][ws.run] is None):
            return '{}/log/{}.{}.fastq_dump.done'.format(config['wrk_folder'], ws.experiment, ws.run)
        else:
            return '{}/log/{}.{}.fastq_gunzip.done'.format(config['wrk_folder'], ws.experiment, ws.run)
    else:
        return '{}/log/{}.{}.trimmomatic.done'.format(config['wrk_folder'], ws.experiment, ws.run)

def get_prefix(ws):
    if (ws.state == 'unfiltered'):
        return '{}/data/{}/{}'.format(config['wrk_folder'], ws.experiment, ws.run)
    else:
        return '{}/data/{}/{}_filtered'.format(config['wrk_folder'], ws.experiment, ws.run)

rule fastqc_all:
    input:
        config['wrk_folder'] + '/log/{experiment}.{run}.unfiltered.fastqc.done',
        config['wrk_folder'] + '/log/{experiment}.{run}.filtered.fastqc.done',

rule fastqc:
    input:      dump_or_gunzip
    output:     touch(config['wrk_folder'] + '/log/{experiment}.{run}.{state}.fastqc.done')
    params:     prefix = get_prefix, 
                layout = lambda ws: config['experiments'][ws.experiment]['layout'],
                outdir =  config['wrk_folder'] + '/qc/{experiment}/{run}_{state}'
    threads:    config['tools']['fastqc']['threads']
    benchmark:  'log/{experiment}.{run}.{state}.fastqc.benchmark'
    log:        outfile = 'log/{experiment}.{run}.{state}.fastqc.out',
                errfile = 'log/{experiment}.{run}.{state}.fastqc.err'
    conda:      config['conda_env_config']
    shell:      '''
                mkdir -p {params.outdir}
                if [[ "{params.layout}" == "PE" ]]; then
                    {config[tools][fastqc][call]}   \
                            -t {threads}            \
                            -o {params.outdir}      \
                            --extract               \
                            -f fastq                \
                            {params.prefix}_1.fastq {params.prefix}_2.fastq \
                            1> {log.outfile} 2> {log.errfile}
                else
                    {config[tools][fastqc][call]}   \
                            -t {threads}            \
                            -o {params.outdir}      \
                            --extract               \
                            -f fastq                \
                            {params.prefix}_1.fastq \
                            1> {log.outfile} 2> {log.errfile}
                fi
                '''
