#!/usr/bin/env python3
"""
Cluster contigs using sklearn.

"""

import sys
import argparse
import logging
from contextlib import redirect_stdout

import pandas as pd
import numpy as np
from sklearn.mixture import GaussianMixture, BayesianGaussianMixture

logger = logging.getLogger(__name__)


def load_data(data_path, length_path, min_length):
    # Load from file.
    data = pd.read_table(args.data_path, index_col='contig_id')
    length = pd.read_table(args.length_path,
                           names=['contig_id', 'length'],
                           index_col='contig_id', squeeze=True)
    # Align indices.
    length = length.loc[data.index]
    # Length filter.
    data = data[length >= min_length]
    length = length.loc[data.index]
    logger.info('Loaded {} contigs with {} bp total.'
                .format(len(length), length.sum()))
    return data, length


def _fit_gmm(data, max_nbins, alpha=None, seed=None):
    logger.info('Fitting a VBGMM model:'
                f' max_nbins={max_nbins}, alpha={alpha}, seed={seed}')
    model = BayesianGaussianMixture(max_nbins, covariance_type='full',
                                    weight_concentration_prior_type='dirichlet_process',
                                    weight_concentration_prior=alpha,
                                    random_state=seed,
                                    max_iter=1000,
                                    verbose=2, verbose_interval=1)
    with redirect_stdout(sys.stderr):
        model.fit(data)
    if model.converged_:
        logger.info('Finished clustering.')
    else:
        logger.warning('Convergence not achieved.'
                       ' Results may not be correct.')
    return model


def _label_gmm(model, data, pthresh=0.0):
    logger.info(f'Calculating posterior probabilities.')
    if pthresh == 0:
        logger.info(f'Assigning bins to all contigs.')
        best_guess = model.predict(data)
        mask = np.ones_like(best_guess).astype(bool)
    elif pthresh > 0:
        logger.info(f'Assigning bins to contigs where p > {pthresh}.')
        probs = model.predict_proba(data)
        best_guess = np.argmax(probs, axis=1)
        # Are any of the probs > pthresh?
        mask = (probs > pthresh).sum(axis=1).astype(bool)
    cluster = pd.Series(best_guess.astype(int),
                        index=data.index, name='cluster')
    return cluster[mask]


# def cluster_gmm(data, max_nbins, alpha=None, seed=None, prob_min=0):
#     model = _fit_gmm(data, max_nbins, alpha, seed)
#     cluster = _label_gmm(model, data, pthresh=prob_min)
#     return cluster
#
#
def cluster_gmm(data, length, max_nbins, frac,
                alpha=None, seed=None, prob_min=0):
    if 0 < frac < 1:
        subdata = data.sample(frac=frac, weights=length, random_state=seed)
        logger.info(f'Subsampling data weighted by contig length: frac={frac}')
        logger.info('Subsampled {} reads with {} bp total.'
                    .format(subdata.shape[0], length.loc[subdata.index].sum()))
    elif frac == 1:
        subdata = data.copy()
    else:
        raise ValueError("*frac* must be in range (0, 1]")
    model = _fit_gmm(subdata, max_nbins, alpha, seed)
    cluster = _label_gmm(model, data, pthresh=prob_min)
    return cluster


def filt_by_total_size(cluster, length, min_bin_size):
    cluster_length = length.groupby(cluster).sum()
    valid_clusters = ((cluster_length >= min_bin_size)
                      [lambda x: x]
                      .index)
    filt_cluster = cluster[cluster.isin(valid_clusters)]
    total_clusters = len(valid_clusters)
    logger.info(f'{total_clusters} clusters with more than {min_bin_size} bp')
    return filt_cluster


# TODO: Institute this.
def rename_clusters(cluster, length):
    cluster_length = length.groupby(cluster).sum()
    rename_by = (-cluster_length).argsort()
    renamed_cluster = cluster.map(rename_by)
    return renamed_cluster 


def summarize_clusters(cluster, length):
    cluster_summary = (length.groupby(cluster).agg(['count', 'sum'])
                             .sort_values('sum', ascending=False))
    cluster_summary.rename(columns={'count': 'n_contigs',
                                    'sum': 'total_length'},
                           inplace=True)
    cluster_summary.rename(index=int, inplace=True)
    logger.info('{} bp in {} valid clusters'
                .format(cluster_summary.total_length.sum(),
                        cluster_summary.shape[0])
                )
    logger.info('Top 10 clusters:\n\n{}'.format(cluster_summary.head(10)))
    return cluster_summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('data_path', metavar='CONCOCT_PCA_TRANSFORM_CSV')
    parser.add_argument('length_path', metavar='CONTIG_LENGTH_TSV')
    parser.add_argument('--min-length', type=int, default=0)
    parser.add_argument('--seed', type=int, default=1)
    parser.add_argument('--verbosity', default='DEBUG')
    parser.add_argument('--summary',
                        help='print summary of clusters to path')

    parser.add_argument('--max-nbins', type=int, default=1000)
    parser.add_argument('--prob-min', type=float, default=0.0,
                        help=('posterior probability minimum to be assigned to'
                              ' a bin'))
    parser.add_argument('--frac', type=float, default=1,
                        help='fraction of data used to fit initial GMM')
    parser.add_argument('--alpha', type=float, default=0,
                        help=('concentration parameter on number of'
                              'componenets; defaults to 0 meaning 1/MAX_NBINS')
                        )
    parser.add_argument('--outfile', '-o', type=argparse.FileType('w'),
                        default=sys.stdout)

    args = parser.parse_args()

    logging.basicConfig(format='%(message)s',
                        level=getattr(logging, args.verbosity.upper()))

    data, length = load_data(args.data_path, args.length_path, args.min_length)

    # Anything <= 0 is nonsensicle for concentration parameter, so we
    # use these as a placeholder for the sklearn default (1 / n_components).
    if args.alpha <= 0:
        args.alpha = None

    clusters = cluster_gmm(data, length,
                           max_nbins=args.max_nbins,
                           alpha=args.alpha,
                           seed=args.seed,
                           frac=args.frac,
                           prob_min=args.prob_min,
                           )

    clusters = rename_clusters(clusters, length)
    cluster_summary = summarize_clusters(clusters, length)

    if args.summary:
        cluster_summary.to_csv(args.summary, sep='\t')
    clusters.to_csv(args.outfile, sep='\t', header=True)
