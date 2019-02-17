
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

rule parse_ncbi_genomes_table:
    output: genome='meta/genome.tsv', replicon='meta/replicon.tsv'
    input: script='scripts/parse_ncbi_genome_tables.py', ncbi='meta/ncbi_genomes.csv'
    shell: '{input.script} {input.ncbi} {output.genome} {output.replicon}'
