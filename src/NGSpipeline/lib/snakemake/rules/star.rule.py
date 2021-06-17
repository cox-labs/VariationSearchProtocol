rule star_index:
    input:      config['ref_folder'] + '/log/{genome}.fa.wget_genome.done',
                config['ref_folder'] + '/log/{genome}.gtf.wget_genome.done'
    output:     touch(config['ref_folder'] + '/log/{genome}.star_index.done')
    params:     input_fasta = config['ref_folder'] + '/{genome}/genome.fa',
                input_gtf = config['ref_folder'] + '/{genome}/genome.gtf',
                genome_dir = config['ref_folder'] + '/{genome}/star_index'
    threads:    config['tools']['star']['index']['threads']
    benchmark:  'log/{genome}.star_index.benchmark'
    log:        outfile = 'log/{genome}.star_index.out',
                errfile = 'log/{genome}.star_index.err'
    conda:      config['conda_env_config']
    shell:      '''
                mkdir -p {params.genome_dir}
                {config[tools][star][index][call]}          \
                    --runThreadN {threads}                  \
                    --genomeDir {params.genome_dir}         \
                    --genomeFastaFiles {params.input_fasta} \
                    --sjdbGTFfile {params.input_gtf}        \
                    1> {log.outfile} 2> {log.errfile}
                '''
################################################################################
#   In order to support loading of large genomes
#       sysctl -w kernel.shmall=8000000
#       sysctl -w kernel.shmmax=32000000000
#       To make the changes permanent, you can change the /etc/sysctl.conf file.
################################################################################

rule star_load:
    input:      config['ref_folder'] + '/log/{genome}.star_index.done'
    output:     touch(config['wrk_folder'] + '/log/{genome}.star_load.done')
    params:     in_memory = lambda ws: '{}/log/{}.star_load.in_memory.{}'.format(config['ref_folder'], ws.genome, hashlib.md5(config['wrk_folder'].encode()).hexdigest()),
                genome_dir = config['ref_folder'] + '/{genome}/star_index/'
    threads:    config['tools']['star']['load']['threads']
    benchmark:  'log/{genome}.star_load.benchmark'
    log:        outfile = 'log/{genome}.star_load.out',
                errfile = 'log/{genome}.star_load.err'
    conda:      config['conda_env_config'] 
    shell:      '''
                touch {params.in_memory}
                {config[tools][star][load][call]}   \
                    --genomeDir {params.genome_dir} \
                    1> {log.outfile} 2> {log.errfile}
                '''

def get_input_fasta(ws):
    if (config['experiments'][ws.experiment]['strategy'] == 'Ribo-Seq'):
        if (config['experiments'][ws.experiment]['runs'][ws.run] is None):
            return '{}/log/{}.{}.fastq_dump.done'.format(config['wrk_folder'], ws.experiment, ws.run)
        else:
            return '{}/log/{}.{}.fastq_gunzip.done'.format(config['wrk_folder'], ws.experiment, ws.run)

    else:
        return '{}/log/{}.{}.trimmomatic.done'.format(config['wrk_folder'], ws.experiment, ws.run)

rule star_align:
    input:      lambda ws: '{}/log/{}.star_load.done'.format(config['wrk_folder'], config['experiments'][ws.experiment]['organism']),
                get_input_fasta
    output:     touch(config['wrk_folder'] + '/log/{experiment}.{run}.star_align.done')
    params:     prefix = config['wrk_folder'] + '/data/{experiment}/{run}',
                layout = lambda ws: config['experiments'][ws.experiment]['layout'],
                genome_dir = lambda ws: '{}/{}/star_index/'.format(config['ref_folder'], config['experiments'][ws.experiment]['organism']),
                strategy = lambda ws: config['experiments'][ws.experiment]['strategy']
    threads:    config['tools']['star']['align']['threads']
    benchmark:  'log/{experiment}.{run}.star_align.benchmark'
    log:        outfile = 'log/{experiment}.{run}.star_align.out',
                errfile = 'log/{experiment}.{run}.star_align.err'
    conda:      config['conda_env_config']
    shell:      '''
                if [[ "{params.layout}" == "PE" ]]; then
                    {config[tools][star][align][call]}                  \
                        --runThreadN {threads}                          \
                        --genomeDir {params.genome_dir}                 \
                        --readFilesIn {params.prefix}_filtered_1.fastq  \
                        {params.prefix}_filtered_2.fastq                \
                        --outFileNamePrefix {params.prefix}.            \
                        1> {log.outfile} 2> {log.errfile}
                else
                    if [[ "{params.strategy}" == "Ribo-Seq" ]]; then
                        {config[tools][star][align][call]}                  \
                            --runThreadN {threads}                          \
                            --genomeDir {params.genome_dir}                 \
                            --readFilesIn {params.prefix}_1.fastq           \
                            --outFileNamePrefix {params.prefix}.            \
                            --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAATCTCGTATG \
                            --alignIntronMin    10                          \
                            --alignMatesGapMax  5000                        \
                            1> {log.outfile} 2> {log.errfile}
                    else
                        {config[tools][star][align][call]}                  \
                            --runThreadN {threads}                          \
                            --genomeDir {params.genome_dir}                 \
                            --readFilesIn {params.prefix}_filtered_1.fastq  \
                            --outFileNamePrefix {params.prefix}.            \
                            1> {log.outfile} 2> {log.errfile}
                    fi
                fi
                mv {params.prefix}.Aligned.out.sam {params.prefix}.sam
                rm {params.prefix}.Log.final.out {params.prefix}.Log.out {params.prefix}.Log.progress.out {params.prefix}.SJ.out.tab
                '''

def get_star_align_dones(ws):
    result = []
    for exp in config['experiments']:
        if (config['experiments'][exp]['organism'] == ws.genome):
            for run in config['experiments'][exp]['runs']:
                result.append('{}/log/{}.{}.star_align.done'.format(config['wrk_folder'], exp, run))
    return result

import hashlib
rule star_remove:
    input:      get_star_align_dones
    output:     touch(config['wrk_folder'] + '/log/{genome}.star_remove.done')
    params:     in_memory = lambda ws: '{}/log/{}.star_load.in_memory.{}'.format(config['ref_folder'], ws.genome, hashlib.md5(config['wrk_folder'].encode()).hexdigest()),
                prefix = lambda ws: '{}/log/{}.star_load.in_memory'.format(config['ref_folder'], ws.genome),
                genome_dir = config['ref_folder'] + '/{genome}/star_index/'
    threads:    config['tools']['star']['remove']['threads']
    benchmark:  'log/{genome}.star_remove.benchmark'
    log:        outfile = 'log/{genome}.star_remove.out',
                errfile = 'log/{genome}.star_remove.err'
    conda:      config['conda_env_config'] 
    shell:      '''
                rm {params.in_memory}
                if [[ `ls {params.prefix}* | wc -l` == "0" ]]; then
                    {config[tools][star][remove][call]} \
                        --genomeDir {params.genome_dir} \
                        1> {log.outfile} 2> {log.errfile}
                fi
                '''
