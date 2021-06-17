import yaml
import os

def con2env(config):
    env_backbone = {'channels':['bioconda', 'defaults', 'conda-forge'], 'dependencies':[]}
    with open(config['conda_env_config'], 'w') as env_stream:
        tools = config['tools']
        env_backbone['dependencies'] = [tool+"="+tools[tool]['version'] for tool in tools]
        yaml.dump(env_backbone, env_stream, default_flow_style=False)


include: 'lib/snakemake/rules/sratools.rule.py'
include: 'lib/snakemake/rules/custom.rule.py'
include: 'lib/snakemake/rules/fastqc.rule.py'
include: 'lib/snakemake/rules/trimmomatic.rule.py'
include: 'lib/snakemake/rules/star.rule.py'
include: 'lib/snakemake/rules/samtools.rule.py'

workdir: config['wrk_folder']
con2env(config) #create env file

def get_fq():
    l = []
    for exp in config['experiments']:
        for run in config['experiments'][exp]['runs']:
            l.append('{}/log/{}.{}.unfiltered.fastqc.done'.format(config['wrk_folder'], exp, run))
            l.append('{}/log/{}.{}.filtered.fastqc.done'.format(config['wrk_folder'], exp, run))
    return l

def get_genomes():
    l = []
    for exp in config['experiments']:
        l.append(config['experiments'][exp]['organism'])
    return ['{}/log/{}.star_remove.done'.format(config['wrk_folder'], org) for org in set(l)]


rule all:
    input:
        ['{}/log/{}.index.done'.format(config['wrk_folder'], exp) for exp in config['experiments']],
        get_fq()

