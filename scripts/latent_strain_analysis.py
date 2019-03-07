#!/usr/bin/env python3

if __name__ == "__main__":
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
    ntotal = nmapping.sum()
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

    np.random.seed(2)

    nu0 = np.random.multinomial(1, np.ones(n_taxa) / n_taxa, size=max_genomes)
    assert not (nu0.sum(0) == 0).any(), 'The random intial value for nu fails to account for all genomes.'
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
        nu_idx = pm.Categorical('nu', p=np.ones(max_genomes) / max_genomes, shape=n_taxa)
        nu = tt.extra_ops.to_one_hot(nu_idx, nb_class=max_genomes).T

        # Latent genome contig content
        # TODO: density_param hyper-prior
        # TODO: Consider making theta continuous (non-negative) instead of
        # discrete.
        theta = pm.Poisson('theta', np.ones((max_genomes, n_contigs)) * density_param,
                           shape=(max_genomes, n_contigs), testval=theta0)

        # Measurement of abundances
        # TODO: Since rabund will sometimes be from 16S data, I need a copy number
        # (and sequencing bias) parameter here.
        expect_rabund = simplex_normalize(pi.dot(nu))
        obs_rabund = pm.Dirichlet('obs_rabund', a=expect_rabund * rabund_precision,
                                  shape=(n_samples, max_genomes), observed=rabund.values)

        # Measurement of depth
        expect_frac_nmapping = simplex_normalize(pi.dot(theta) * nlength.values)
        # TODO: Consider using an over-dispersed multinomial so this isn't
        # _too_ influential.
        obs_nmapping = pm.Multinomial('obs_nmapping', n=ntotal.values, p=expect_frac_nmapping,
                                      shape=(n_contigs, n_samples), observed=nmapping.values)
