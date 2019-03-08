#!/usr/bin/env python3
import pandas as pd
import pymc3 as pm
import matplotlib.pyplot as plt
import numpy as np
import theano.tensor as tt

if __name__ == "__main__":
    # TODO: Take these as input params
    cvrg_path = 'data/sim/escherichia_congenic.n1e6.z16.s00.a.backmap.cvrg.tsv'  # sys.argv[1]
    nlength_path = 'data/sim/escherichia_congenic.n1e6.z16.s00.a.nlength.tsv'  # sys.argv[2]
    abund_path = 'data/sim/escherichia_congenic.n1e6.z16.s00.a.abund.tsv'  # sys.argv[3]
    cvrg  = pd.read_table(cvrg_path,
                          names=['sample_id', 'contig_id', 'coverage'],
                          index_col=['sample_id', 'contig_id'], squeeze=True
                         ).unstack('sample_id', fill_value=0)
    nlength = pd.read_table(nlength_path,
                            names=['contig_id', 'nlength'], index_col=['contig_id'], squeeze=True)
    nmapping = cvrg.multiply(nlength, axis=0).round().astype(int).T
    ntotal = nmapping.sum(1)
    abund = pd.read_table(abund_path,
                          names=['sample_id', 'genome_id', 'abundance'],
                          index_col=['sample_id', 'genome_id'], squeeze=True
                         ).unstack('genome_id', fill_value=0)

    rabund = abund.divide(abund.sum(1), axis=0)


    n_taxa = len(rabund.columns)
    n_samples = len(cvrg.columns)
    n_contigs = len(cvrg.index)

    # TODO: Take these as input params
    max_genomes = 10
    rabund_precision = 10000
    diversity_param = 1  # Regularization of numbers of latent components.
    density_param = 0.5  # What fraction of the contigs are expected to be in each genome.
    rabund_noise_param = 1e-6  # This is used to stabilize the calculation, since the llk
                               # is undefined if one of the measured taxa is not related to
                               # to any underlying genome.

    np.random.seed(10)

    #nu0 = np.random.multinomial(1, np.ones(n_taxa) / n_taxa, size=max_genomes)
    nu_idx0 = np.random.choice(range(n_taxa), size=max_genomes)
    nu0 = tt.extra_ops.to_one_hot(nu_idx0, nb_class=n_taxa).eval()
    # assert not (nu0.sum(0) == 0).any()
    theta0 = np.random.multinomial(1, np.ones(max_genomes) / max_genomes, size=n_contigs).T

    with pm.Model() as model0:
        # Latent genome relative abundances
        # TODO: diversity_param hyper-prior
        pi = pm.Dirichlet('pi', a=np.ones(max_genomes) * diversity_param,
                          shape=(n_samples, max_genomes))

        # Latent genome taxonomic identity
        # TODO: Are all latent genomes reflected in the taxonomic measurements?
        # nu = pm.Multinomial('nu', n=0, p=np.ones((max_genomes, n_taxa)) / max_genomes,
        #                     shape=(max_genomes, n_taxa), testval=nu0)
        nu_idx = pm.Categorical('nu_idx', p=np.ones(n_taxa) / n_taxa, shape=max_genomes, testval=nu_idx0)
        nu = tt.extra_ops.to_one_hot(nu_idx, nb_class=n_taxa)


        # Latent genome contig content
        # TODO: density_param hyper-prior
        # TODO: Consider making theta continuous (non-negative) instead of
        # discrete.
        theta = pm.Poisson('theta', np.ones((max_genomes, n_contigs)) * density_param,
                           shape=(max_genomes, n_contigs), testval=theta0)

        # Measurement of abundances
        # TODO: Since rabund will sometimes be from 16S data, I need a copy number
        # (and sequencing bias) parameter here.
        expect_rabund = simplex_normalize(pi.dot(nu) + rabund_noise_param)
        obs_rabund = pm.Dirichlet('obs_rabund', a=expect_rabund * rabund_precision,
                                  shape=(n_samples, max_genomes), observed=rabund.values)

        # Measurement of depth
        expect_frac_nmapping = simplex_normalize(pi.dot(theta) * nlength.values)
        # TODO: Consider using an over-dispersed multinomial so this isn't
        # _too_ influential.
        obs_nmapping = pm.Multinomial('obs_nmapping', n=ntotal.values, p=expect_frac_nmapping,
                                      shape=(n_contigs, n_samples), observed=nmapping.values)

    with model0:
        trace0 = pm.sample()
