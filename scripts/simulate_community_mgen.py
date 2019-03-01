#!/usr/bin/env python3

import sys
import os.path
import pandas as pd
import numpy as np
from tqdm import tqdm
import iss
from iss.generator import simulate_read
from iss.error_models.kde import KDErrorModel
from iss.util import convert_n_reads
from Bio.SeqIO import index as seq_file_index
from Bio.SeqIO import write as write_seq

PLASMID_MULTI = 10
ERR_MODEL_NAME = 'HiSeq'

if __name__ == '__main__':
    community_path = sys.argv[1]
    community_id = sys.argv[2]
    replicon_path = sys.argv[3]
    genome_dir = sys.argv[4]
    nreads = int(float(sys.argv[5]))
    seed = int(sys.argv[6])
    abund_path = sys.argv[7]
    r1_path = sys.argv[8]
    r2_path = sys.argv[9]

    print('Loading input data.', file=sys.stderr)
    community = (pd.read_table(community_path)
                     [lambda x: x.community_id == community_id])
    community.set_index('genome_id', inplace=True)
    assert community.index.is_unique, \
            f'Some genomes are repeated in the community table.'
    replicon = (pd.read_table(replicon_path)
                    [lambda x: x.genome_id.isin(community.index)])
    replicon.set_index('replicon_id', inplace=True)
    assert replicon.index.is_unique, \
            'Some replicons are found multiple times in the replicon table.'
    assert set(replicon.genome_id) == set(community.index), \
            'Not all community genomes are found in the replicon table.'

    print('Simulating read counts.', file=sys.stderr)
    np.random.seed(seed)
    community['sim_abund'] = np.random.lognormal(mean=community['mean_log_abund'],
                                                 sigma=community['stdev_log_abund'])
    community['sim_abund'].to_csv(abund_path, sep='\t')
    replicon['sim_copies'] = (replicon['genome_id'].map(community['sim_abund']) *
                              replicon['replicon_type']
                                      .map({'chromosome': 1,
                                            'plasmid': PLASMID_MULTI})
                                      .fillna(1))
    replicon['sim_nucleotides'] = replicon['sim_copies'] * replicon['size']
    replicon['sim_pread'] = (replicon['sim_nucleotides'] /
                             replicon['sim_nucleotides'].sum())
    replicon['sim_nreads'] = np.random.multinomial(nreads, replicon['sim_pread'])
    print(replicon[['genbank_id', 'size', 'sim_nreads']], file=sys.stderr)

    print('Simulating metagenome reads.', file=sys.stderr)
    model_path = os.path.join(iss.__path__[0], 'profiles', ERR_MODEL_NAME)
    err_model = KDErrorModel(model_path)
    with open(r1_path, 'w') as handle1, open(r2_path, 'w') as handle2:
        reads_so_far = 0
        for genome_id, g in replicon.groupby('genome_id'):
            print(f'Simulating reads for genome {genome_id}', file=sys.stderr)
            genome_path = os.path.join(genome_dir, genome_id + '.fn')
            seqs = seq_file_index(genome_path, 'fasta')
            for replicon_id, r in g.iterrows():
                print(f'Simulating reads for replicon {replicon_id}', file=sys.stderr)
                record = seqs[r.genbank_id]
                for i in tqdm(range(r.sim_nreads), initial=reads_so_far, total=nreads):
                    paired_reads = simulate_read(record, err_model, i)
                    write_seq(paired_reads[0], handle1, 'fastq-sanger')
                    write_seq(paired_reads[1], handle2, 'fastq-sanger')
                reads_so_far += r.sim_nreads
