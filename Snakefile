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
    output: "data/{library}.m.k{klen}.jf"
    input: r1="data/{library}.m.r1.fq.gz", r2="data/{library}.m.r2.fq.gz"
    threads: 2
    shell:
        """
        jellyfish count -t {threads} -m {wildcards.klen} -s 10M -C -o {output} <(zcat {input.r1}) <(zcat {input.r2})
        """

rule histo_counts:
    output: "data/{stem}.histo"
    input: "data/{stem}.jf"
    shell: "jellyfish histo -o {output} {input}"

rule dump_counts:
    output: "data/{stem}.tally"
    input: "data/{stem}.jf"
    shell:
        """
        tmp=`mktemp`
        jellyfish dump -c {input} > $tmp
        sort $tmp > {output}
        """

rule test_sparse_svd:
    output: "build/{asmbl_group}.a.k{klen}.test_sparse_svd.tsv"
    input: script="scripts/test_sparse_svd.py",
        tallies=lambda w: [f'data/{library}.m.k{w.klen}.tally' for library in config['asmbl_group'][w.asmbl_group]]
    params: lib_list=lambda w: [f'{library}:data/{library}.m.k{w.klen}.tally' for library in config['asmbl_group'][w.asmbl_group]]
    shell:
        """
        {input.script} {params.lib_list} > {output}
        """

rule combine_counts:
    output: "data/{asmbl_group}.a.k{klen}.counts"
    input: lambda w: [f'data/{library}.m.k{w.klen}.tally' for library in config['asmbl_group'][w.asmbl_group]]
    params: libs=lambda w: [library for library in config['asmbl_group'][w.asmbl_group]]
    shell:
        """
        rm -rf {output}
        for lib in {params.libs}
        do
            awk -v OFS='\t' -v lib=$lib '{{print lib,$1,$2}}' data/$lib.m.k{wildcards.klen}.tally >> {output}
        done
        """



