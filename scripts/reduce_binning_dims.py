#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

READ_LENGTH = 100
EXPLAIN_THRESH = 0.9
FLOAT_FORMAT = '%.6g'

if __name__ == '__main__':
    length_path = sys.argv[1]
    kmer_path = sys.argv[2]
    cvrg_path = sys.argv[3]
    min_contig_length = int(sys.argv[4])

    print('Loading data', file=sys.stderr)
    length = pd.read_table(length_path, names=['contig_id', 'length'],
                           index_col='contig_id',
                           squeeze=True)
    kmer = pd.read_table(kmer_path, names=['contig_id', 'kmer', 'tally'],
                         index_col=['contig_id', 'kmer'],
                         squeeze=True)
    cvrg = pd.read_table(cvrg_path,
                         names=['library_id', 'contig_id', 'coverage'],
                         index_col=['contig_id', 'library_id'],
                         squeeze=True)

    print('Reshaping data', file=sys.stderr)
    kmer = kmer.unstack('kmer', fill_value=0)
    cvrg = cvrg.unstack('library_id', fill_value=0)
    cvrg = cvrg.reindex(length.index).fillna(0)  # Some contigs without coverage

    print('Transforming data', file=sys.stderr)
    kmer_ps = kmer + 1  # Pseudocount
    cvrg_ps = cvrg.add(READ_LENGTH / length, axis='index')  # Pseudocount
    kmer_comp = kmer_ps.div(length, axis='index')
    # TODO: Consider fixing this so that coverage is normalized correctly across
    # larger/smaller contig lengths.
    # Simplex transform over libraries and then over contigs
    cvrg_norm1 = (cvrg_ps / cvrg_ps.sum())
    cvrg_norm2 = cvrg_norm1.div(cvrg_norm1.sum(1), axis='index')
    total_cvrg = cvrg_ps.sum(axis='columns')

    table = pd.concat([kmer_comp, cvrg_norm2], sort=True, axis='columns')
    table['total_cvrg'] = total_cvrg
    norm_log_table = np.log(table).apply(lambda x: (x - x.mean()) / x.std())

    print('Dimensional reduction', file=sys.stderr)
    # Only fit on long-enough contigs
    train_data = norm_log_table[length > min_contig_length]  # FIXME: This filter throws a warning
    pca = PCA().fit(train_data)
    # Transform all contigs
    num_components = (pca.explained_variance_ratio_.cumsum() < EXPLAIN_THRESH).sum() + 1
    out = pd.DataFrame(pca.transform(norm_log_table)[:,:num_components],
                       index=length.index,
                       columns=['PC{}'.format(i + 1) for i in range(num_components)])
    out.to_csv(sys.stdout, sep='\t', float_format=FLOAT_FORMAT)
