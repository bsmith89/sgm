
MAX_THREADS = 24

include: 'snake/real_data.snake'
include: 'snake/kmer.snake'
include: 'snake/simulate_mgen.snake'

rule all:
    output: []
