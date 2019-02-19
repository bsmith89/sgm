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
        echo "Please activate your environment and then run `pip install -r requirements.txt` or analogous."
        '''

rule start_jupyter:
    shell: 'jupyter notebook --config=nb/jupyter_notebook_config.py --notebook-dir=nb/'

rule download_genbank_bac_table:
    output: 'raw/bacterial_assemblies.tsv'
    params:
        url='ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt'
    resources:
        net_requests_per_second=1
    shell:
        curl_recipe

rule parse_genbank_genomes_table:
    output: genome='data/genome.tsv', replicon='data/replicon.tsv'
    input:
        script='scripts/parse_genbank_genome_table.py',
        table='raw/bacterial_assemblies.tsv',
        assmbl_rprt='raw/ref/genbank_genome_assembly_reports/'
    shell:
        '''
        {input.script} {input.table} {input.assmbl_rprt} {output.genome} {output.replicon}
        '''
