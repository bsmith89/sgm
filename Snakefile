from lib.snake import curl_recipe

MAX_THREADS = 24

include: 'snake/real_data.snake'
include: 'snake/kmer.snake'
include: 'snake/simulate_mgen.snake'

rule all:
    output: []

rule initialize_project:
    shell:
        '''
        git config --local filter.dropoutput_ipynb.clean scripts/ipynb_output_filter.py
        git config --local filter.dropoutput_ipynb.smudge cat
        git config --local diff.daff-csv.command "daff.py diff --git"
        git config --local merge.daff-csv.name "daff.py tabular merge"
        git config --local merge.daff-csv.driver "daff.py merge --output %A %O %A %B"
        echo "Please activate your environment and then run `pip install -r requirements.txt` or analagous."
        '''

rule start_jupyter:
    shell: 'jupyter notebook --config=nb/jupyter_notebook_config.py --notebook-dir=nb/'

rule download_genbank_bac_table:
    output: 'raw/bacterial_assemblies.tsv'
    params:
        url='ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'
    shell:
        curl_recipe

rule parse_genbank_genomes_table:
    output: genome='data/genome.tsv'
    input: script='scripts/parse_genbank_genome_table.py', table='raw/bacterial_assemblies.tsv'
    shell: '{input.script} {input.table} > {output.genome}'
