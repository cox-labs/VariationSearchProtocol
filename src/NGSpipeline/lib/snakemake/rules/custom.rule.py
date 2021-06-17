rule wget_genome:
    input:
    output:     touch(config['ref_folder'] + '/log/{genome}.{type}.wget_genome.done')
    params:     ftp_address = lambda ws: config['resources']['genomes'][ws.genome][ws.type],
                genome_file = config['ref_folder'] + '/{genome}/genome.{type}',
                folder = config['ref_folder'] + '/{genome}/'
    threads:    8
    benchmark:  config['ref_folder'] + '/log/{genome}.{type}.wget_genome.benchmark'
    log:        outfile = config['ref_folder'] + '/log/{genome}.{type}.wget_genome.out',
                errfile = config['ref_folder'] + '/log/{genome}.{type}.wget_genome.err'
    shell:      '''
                mkdir -p {params.folder}
                wget {params.ftp_address} -O {params.genome_file}.gz 1> {log.outfile} 2> {log.errfile}
                gunzip {params.genome_file}.gz 1> {log.outfile} 2> {log.errfile}
                '''

rule get_adapters:
    input:
    output:     touch(config['ref_folder'] + '/log/get_adapters.done')
    params:     adapter_in = config['lib_folder'] + '/lib/adapters/adapters.fa',
                folder_out = config['ref_folder'] + '/'
    threads:    1
    benchmark:  config['ref_folder'] + '/log/get_adapters.benchmark'
    log:        outfile = config['ref_folder'] + '/log/get_adapters.out',
                errfile = config['ref_folder'] + '/log/get_adapters.err'
    shell:      '''
                mkdir -p {params.folder_out}
                cp {params.adapter_in} {params.folder_out}. 1> {log.outfile} 2> {log.errfile}
                '''

rule gunzip_to_wrk_dir:
    input:      
    output:     touch(config['wrk_folder'] + '/log/{experiment}.{run}.fastq_gunzip.done')
    params:     prefix = config['wrk_folder'] + '/data/{experiment}/{run}',
                files = lambda ws: config['experiments'][ws.experiment]['runs'][ws.run],
                layout = lambda ws: config['experiments'][ws.experiment]['layout'],
                folder = config['wrk_folder'] + '/data/{experiment}/'
    threads:    8
    benchmark:  'log/{experiment}.{run}.gunzip_to_wrk_dir.benchmark'
    log:        outfile = 'log/{experiment}.{run}.gunzip_to_wrk_dir.out',
                errfile = 'log/{experiment}.{run}.gunzip_to_wrk_dir.err'
    run:
                shell('mkdir -p {params.folder}')
                for i in range(len(params.files)):
                    shell('gunzip -c ' + params.files[i] + ' 1> {params.prefix}_' + str(i+1) + '.fastq 2> {log.errfile}')
