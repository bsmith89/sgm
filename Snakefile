import pandas as pd
from lib.snake import alias_recipe, alias_fmt, curl_recipe, curl_unzip_recipe

# Configure the pipeline
config_file = 'config.yaml'
configfile: config_file
# Metadata specific configurations
_library = pd.read_table(config['_meta_library'], index_col='library_id')
_asmbl_group = pd.read_table(config['_meta_asmbl_group'])
config['library'] = {}
for library_id, row in _library.iterrows():
    config['library'][library_id] = {}
    config['library'][library_id]['r1'] = row['file_r1']
    config['library'][library_id]['r2'] = row['file_r2']
config['asmbl_group'] = {}
for group, d in _asmbl_group.groupby('asmbl_group'):
    config['asmbl_group'][group] = list(d['library_id'])


rule all:
    output: []

rule alias_raw_read_r1:
    output: 'data/{library}.m.r1.fq.gz',
    input: lambda wildcards: 'raw/mgen/{}'.format(config['library'][wildcards.library]['r1'])
    shell: alias_recipe

rule alias_raw_read_r2:
    output: 'data/{library}.m.r2.fq.gz',
    input: lambda wildcards: 'raw/mgen/{}'.format(config['library'][wildcards.library]['r2'])
    shell: alias_recipe

rule tally_library:
    output:
        temp("data/{library}.m.k{klen}.kmc_suf"),
        temp("data/{library}.m.k{klen}.kmc_pre"),
    input: r1="data/{library}.m.r1.fq.gz", r2="data/{library}.m.r2.fq.gz"
    threads: 4
    resources:
        mem_gb=10
    shell:
        """
        tmpdir=`mktemp -d`
        tmp=`mktemp`
        echo {input.r1} > $tmp
        echo {input.r2} >> $tmp
        kmc -k{wildcards.klen} -m{resources.mem_gb} -t{threads} -ci1 -cs1000000000 @$tmp data/{wildcards.library}.m.k{wildcards.klen} $tmpdir
        """

# rule histo_counts:
#     output: "data/{stem}.histo"
#     input: "data/{stem}.jf"
#     shell: "jellyfish histo -o {output} {input}"
#
rule dump_counts:
    output: temp("data/{stem}.tally")
    input: "data/{stem}.kmc_suf", "data/{stem}.kmc_pre"
    shell:
        """
        kmc_tools transform data/{wildcards.stem} dump -s {output}
        """

rule merge_counts:
    output: "data/{asmbl_group}.a.k{klen}.counts.tsv"
    input: script="scripts/merge_kmer_counts.py",
        tallies=lambda w: [f'data/{library}.m.k{w.klen}.tally' for library in config['asmbl_group'][w.asmbl_group]]
    params:
        lib_list=lambda w: [f'{library}:data/{library}.m.k{w.klen}.tally' for library in config['asmbl_group'][w.asmbl_group]],
        min_libs = 3,
        min_counts = 5,
    shell:
        """
        {input.script} {params.min_libs} {params.min_counts} {params.lib_list} > {output}
        """

rule simulate_reads:
    output: "sim/test.fq.gz"
    input: "meta/camisim.config.ini"
    conda: "camisim"
    shell:
        '''
        metagenomesimulation.py {input}
        '''
