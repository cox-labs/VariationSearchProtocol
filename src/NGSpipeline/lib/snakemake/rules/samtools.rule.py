aligner_dict = { 
    'WGS' : 'bwa_mem',
    'RNA-Seq' : 'star_align',
    'Ribo-Seq' : 'star_align'
}


def get_input_samtools_sam2bam(ws):
    return '{}/log/{}.{}.{}.done'.format(config['wrk_folder'], ws.experiment, ws.run, aligner_dict[config['experiments'][ws.experiment]['strategy']])

rule samtools_sam2bam:
    input:      get_input_samtools_sam2bam 
    output:	touch(config['wrk_folder'] + '/log/{experiment}.{run}.samtools_sam2bam.done')
    params:     input_sam = config['wrk_folder'] + '/data/{experiment}/{run}.sam',
                output_bam = config['wrk_folder'] + '/data/{experiment}/{run}.unsorted.bam'
    threads:    config['tools']['samtools']['view']['sam2bam']['threads']
    benchmark:  'log/{experiment}.{run}.samtools_sam2bam.benchmark'
    log:        outfile = 'log/{experiment}.{run}.samtools_sam2bam.out',
                errfile = 'log/{experiment}.{run}.samtools_sam2bam.err'
    conda:      config['conda_env_config']
    shell:      '''
                {config[tools][samtools][view][sam2bam][call]}  \
                    -@ {threads}                                \
                    {params.input_sam}                          \
                    -o {params.output_bam}                      \
                    1> {log.outfile} 2> {log.errfile}   
                rm {params.input_sam}
                '''
rule samtools_sort_by_name:
    input:      config['wrk_folder'] + '/log/{experiment}.{run}.samtools_sam2bam.done'
    output:     touch(config['wrk_folder'] + '/log/{experiment}.{run}.samtools_sort_by_name.done')
    params:     param = config['tools']['samtools']['sort']['by_name'],
                input_bam = config['wrk_folder'] + '/data/{experiment}/{run}.unsorted.bam',
                output_bam = config['wrk_folder'] + '/data/{experiment}/{run}.sorted_by_name.bam'
    threads:    config['tools']['samtools']['sort']['threads']
    benchmark:  'log/{experiment}.{run}.samtools_sort_by_name.benchmark'
    log:        outfile = 'log/{experiment}.{run}.samtools_sort_by_name.out',
                errfile = 'log/{experiment}.{run}.samtools_sort_by_name.err'
    conda:      config['conda_env_config']
    shell:      '''
                {config[tools][samtools][sort][call]}   \
                    -@ {threads}                        \
                    {params.param}                      \
                    {params.input_bam}                  \
                    -o {params.output_bam}              \
                    1> {log.outfile} 2> {log.errfile}
                rm {params.input_bam}
                '''

rule samtools_fixmate:
    input:      config['wrk_folder'] + '/log/{experiment}.{run}.samtools_sort_by_name.done'
    output:     touch(config['wrk_folder'] + '/log/{experiment}.{run}.fixmate.done')
    params:     input_bam = config['wrk_folder'] + '/data/{experiment}/{run}.sorted_by_name.bam',
                output_bam = config['wrk_folder'] + '/data/{experiment}/{run}.fixmate.bam'
    threads:    config['tools']['samtools']['fixmate']['threads']
    benchmark:  'log/{experiment}.{run}.samtools_fixmate.benchmark'
    log:        outfile = 'log/{experiment}.{run}.samtools_fixmate.out',
                errfile = 'log/{experiment}.{run}.samtools_fixmate.err'
    conda:      config['conda_env_config']
    shell:      '''
                {config[tools][samtools][fixmate][call]}    \
                    {params.input_bam}                      \
                    {params.output_bam}                     \
                    1> {log.outfile} 2> {log.errfile}
                rm {params.input_bam}
                '''

def get_samtools_sort_by_position_input(ws):
    if (config['experiments'][ws.experiment]['layout'] == 'PE'):
        return ['{}/log/{}.{}.fixmate.done'.format(config['wrk_folder'], ws.experiment, i) for i in config['experiments'][ws.experiment]['runs'].keys()]
    else:
        return ['{}/log/{}.{}.samtools_sam2bam.done'.format(config['wrk_folder'], ws.experiment, i) for i in config['experiments'][ws.experiment]['runs'].keys()]

def get_samtools_sort_by_position_input_bam(ws):
    if (config['experiments'][ws.experiment]['layout'] == 'PE'):
        return ['{}/data/{}/{}.fixmate.bam'.format(config['wrk_folder'], ws.experiment, i) for i in config['experiments'][ws.experiment]['runs'].keys()]
    else:
        return ['{}/data/{}/{}.unsorted.bam'.format(config['wrk_folder'], ws.experiment, i) for i in config['experiments'][ws.experiment]['runs'].keys()]

rule samtools_sort_by_position:
    input:      get_samtools_sort_by_position_input
    output:     touch(config['wrk_folder'] + '/log/{experiment}.{run}.sorted_by_position.done')
    params:     input_bam = get_samtools_sort_by_position_input_bam,
                output_bam = config['wrk_folder'] + '/data/{experiment}/{run}.sorted_by_position.bam'
    threads:    config['tools']['samtools']['sort']['threads']
    benchmark:  'log/{experiment}.{run}.samtools_sort_by_position.benchmark'
    log:        outfile = 'log/{experiment}.{run}.samtools_sort_by_position.out',
                errfile = 'log/{experiment}.{run}.samtools_sort_by_position.err'
    conda:      config['conda_env_config']
    shell:      '''
                {config[tools][samtools][sort][call]}   \
                    -@ {threads}                        \
                    {params.input_bam}                  \
                    -o {params.output_bam}              \
                    1> {log.outfile} 2> {log.errfile}
                rm {params.input_bam}
                '''

rule samtools_markdup:
    input:      config['wrk_folder'] + '/log/{experiment}.{run}.sorted_by_position.done'
    output:     touch(config['wrk_folder'] + '/log/{experiment}.{run}.markdup.done')
    params:     input_bam = config['wrk_folder'] + '/data/{experiment}/{run}.sorted_by_position.bam',
                output_bam = config['wrk_folder'] + '/data/{experiment}/{run}.markdup.bam'
    threads:    config['tools']['samtools']['markdup']['threads']
    benchmark:  'log/{experiment}.{run}.samtools_markdup.benchmark'
    log:        outfile = 'log/{experiment}.{run}.samtools_markdup.out',
                errfile = 'log/{experiment}.{run}.samtools_markdup.err'
    conda:      config['conda_env_config']
    shell:      '''
                {config[tools][samtools][markdup][call]}    \
                    {params.input_bam}                      \
                    {params.output_bam}                     \
                    1> {log.outfile} 2> {log.errfile}
                rm {params.input_bam}
                '''

def get_samtools_merge_input(ws):
    if (config['experiments'][ws.experiment]['layout'] == 'PE'):
        return ['{}/log/{}.{}.markdup.done'.format(config['wrk_folder'], ws.experiment, i) for i in config['experiments'][ws.experiment]['runs'].keys()]
    else:
        return ['{}/log/{}.{}.sorted_by_position.done'.format(config['wrk_folder'], ws.experiment, i) for i in config['experiments'][ws.experiment]['runs'].keys()]

def get_samtools_merge_input_bam(ws):
    if (config['experiments'][ws.experiment]['layout'] == 'PE'):
        return ['{}/data/{}/{}.markdup.bam'.format(config['wrk_folder'], ws.experiment, i) for i in config['experiments'][ws.experiment]['runs'].keys()]
    else:
        return ['{}/data/{}/{}.sorted_by_position.bam'.format(config['wrk_folder'], ws.experiment, i) for i in config['experiments'][ws.experiment]['runs'].keys()]

rule samtools_merge:
    input:      get_samtools_merge_input
    output:     touch(config['wrk_folder'] + '/log/{experiment}.merge.done')
    params:     input_bam = get_samtools_merge_input_bam, 
                output_bam = config['wrk_folder'] + '/data/{experiment}.bam'
    threads:    config['tools']['samtools']['merge']['threads']
    benchmark:  'log/{experiment}.merge.benchmark'
    log:        outfile = 'log/{experiment}.merge.out',
                errfile = 'log/{experiment}.merge.err'
    conda:      config['conda_env_config']
    shell:      '''
                if [[ `echo {params.input_bam} | wc -w` == "1" ]]; then
                    mv {params.input_bam} {params.output_bam}
                else
                    {config[tools][samtools][merge][call]}  \
                        -@ {threads}                        \
                        {params.output_bam}                 \
                        {params.input_bam}                  \
                        1> {log.outfile} 2> {log.errfile}
                    rm {params.input_bam}
                fi
                '''

rule samtools_index:
    input:      config['wrk_folder'] + '/log/{experiment}.merge.done'
    output:     touch(config['wrk_folder'] + '/log/{experiment}.index.done')
    params:     input_bam = config['wrk_folder'] + '/data/{experiment}.bam',
                output_bai = config['wrk_folder'] + '/data/{experiment}.bai'
    threads:    config['tools']['samtools']['index']['threads']
    benchmark:  'log/{experiment}.samtools_index.benchmark'
    log:        outfile = 'log/{experiment}.samtools_index.out',
                errfile = 'log/{experiment}.samtools_index.err'
    conda:      config['conda_env_config']
    shell:      '''
                {config[tools][samtools][index][call]}  \
                    {params.input_bam}                  \
                    {params.output_bai}                 \
                    1> {log.outfile} 2> {log.errfile}
                '''

rule samtools_faidx:
    input:      config['ref_folder'] + '/log/{genome}.fa.wget_genome.done'
    output:     touch(config['ref_folder'] + '/log/{genome}.samtools_faidx.done')
    params:     input_fa = config['ref_folder'] + '/{genome}/genome.fa'
    threads:    config['tools']['samtools']['faidx']['threads']
    benchmark:  config['ref_folder'] + '/log/{genome}.samtools_faidx.benchmark'
    log:        outfile = config['ref_folder'] + '/log/{genome}.samtools_faidx.out',
                errfile = config['ref_folder'] + '/log/{genome}.samtools_faidx.err'
    conda:      config['conda_env_config']
    shell:      '''
                {config[tools][samtools][faidx][call]}  \
                    {params.input_fa}                   \
                    1> {log.outfile} 2> {log.errfile}
                '''

rule bcftools_mpileup:
    input:      lambda ws: '{}/log/{}.fa.wget_genome.done'.format(config['ref_folder'], config['experiments'][ws.experiment]['organism']),
                lambda ws: '{}/log/{}.samtools_faidx.done'.format(config['ref_folder'], config['experiments'][ws.experiment]['organism']),
                config['wrk_folder'] + '/log/{experiment}.merge.done'
    output:     touch(config['wrk_folder'] + '/log/{experiment}.mpileup.done')
    params:     input_fa = lambda ws: '{}/{}/genome.fa'.format(config['ref_folder'], config['experiments'][ws.experiment]['organism']),
                input_bam = config['wrk_folder'] + '/data/{experiment}.bam',
                output_vcf = config['wrk_folder'] + '/data/{experiment}.mpileup.vcf'
    threads:    config['tools']['bcftools']['mpileup']['threads']
    benchmark:  'log/{experiment}.bcftools_mpileup.benchmark'
    log:        errfile_mpileup = 'log/{experiment}.bcftools_mpileup.err',
                outfile_mpileup = 'log/{experiment}.bcftools_mpileup.out'
    conda:      config['conda_env_config']
    shell:      '''
                {config[tools][bcftools][mpileup][call]}    \
                    -f {params.input_fa}                    \
                    {params.input_bam}                      \
                    2> {log.errfile_mpileup} |              \
                    {config[tools][bcftools][call][call]}   \
                    -o {params.output_vcf}                  \
                    1> {log.outfile_mpileup} 2> {log.errfile_mpileup}
                '''
