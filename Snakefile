
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
        echo "Please activate your environment and then run `pip install -r requirements.txt` or analagous."
        '''

rule start_jupyter:
    shell: 'jupyter notebook --config=nb/jupyter_notebook_config.py --notebook-dir=nb/'

