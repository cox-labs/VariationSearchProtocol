def dump_or_gunzip(ws):
    if (config['experiments'][ws.experiment]['runs'][ws.run] is None):
        return '{}/log/{}.{}.fastq_dump.done'.format(config['wrk_folder'], ws.experiment, ws.run)
    else:
        return '{}/log/{}.{}.fastq_gunzip.done'.format(config['wrk_folder'], ws.experiment, ws.run)

rule trimmomatic:
    input:      dump_or_gunzip,
                adapter_do = config['ref_folder'] + '/log/get_adapters.done'
    output:     touch(config['wrk_folder'] + '/log/{experiment}.{run}.trimmomatic.done')
    params:     adapter = config['ref_folder'] + '/adapters.fa',
                prefix = config['wrk_folder'] + '/data/{experiment}/{run}',
                layout = lambda ws: config['experiments'][ws.experiment]['layout'], 
                phred = config['tools']['trimmomatic']['phred'],
                seed_mismatches = config['tools']['trimmomatic']['illumina_clip']['seed_mismatches'],
                palindrome_clip_threshold = config['tools']['trimmomatic']['illumina_clip']['palindrome_clip_threshold'],
                simple_clip_threshold = config['tools']['trimmomatic']['illumina_clip']['simple_clip_threshold'],
                leading = config['tools']['trimmomatic']['leading'],
                trailing = config['tools']['trimmomatic']['trailing'],
                window_size = config['tools']['trimmomatic']['slidingwindow']['window_size'],
                required_quality =config['tools']['trimmomatic']['slidingwindow']['required_quality'],
                min_length = config['tools']['trimmomatic']['minlen']
    threads:    config['tools']['trimmomatic']['threads']
    benchmark:  'log/{experiment}.{run}.trimmomatic.benchmark'
    log:        outfile = 'log/{experiment}.{run}.trimmomatic.out',
                errfile = 'log/{experiment}.{run}.trimmomatic.err'
    conda:      config['conda_env_config']
    shell:      '''
                if [[ "{params.layout}" == "PE" ]]; then
                    {config[tools][trimmomatic][call]}      \
                        {params.layout}                     \
                        {params.phred}                      \
                        -threads {threads}                  \
                        {params.prefix}_1.fastq {params.prefix}_2.fastq \
                        {params.prefix}_filtered_1.fastq {params.prefix}_unfiltered_1.fastq \
                        {params.prefix}_filtered_2.fastq {params.prefix}_unfiltered_2.fastq \
                        ILLUMINACLIP:{params.adapter}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} \
                        LEADING:{params.leading} TRAILING:{params.trailing} \
                        SLIDINGWINDOW:{params.window_size}:{params.required_quality} \
                        MINLEN:{params.min_length}          \
                        1> {log.outfile} 2> {log.errfile}
                    rm {params.prefix}_unfiltered_1.fastq {params.prefix}_unfiltered_2.fastq
                else
                    {config[tools][trimmomatic][call]}      \
                        {params.layout}                     \
                        {params.phred}                      \
                        -threads {threads}                  \
                        {params.prefix}_1.fastq             \
                        {params.prefix}_filtered_1.fastq    \
                        ILLUMINACLIP:{params.adapter}:{params.seed_mismatches}:{params.palindrome_clip_threshold}:{params.simple_clip_threshold} \
                        LEADING:{params.leading} TRAILING:{params.trailing} \
                        SLIDINGWINDOW:{params.window_size}:{params.required_quality} \
                        MINLEN:{params.min_length}          \
                        1> {log.outfile} 2> {log.errfile}
                fi
                '''
