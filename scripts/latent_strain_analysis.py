#!/usr/bin/env python3
import pandas as pd
import pymc3 as pm
import matplotlib.pyplot as plt
import numpy as np
import theano.tensor as tt
import sys

tt_simplex_normalize = lambda x: x / x.sum(1).reshape((x.shape[0], 1))

info = lambda s, *args: print(s % args, file=sys.stderr)

if __name__ == "__main__":
    seq_cvrg_path = sys.argv[1]
    nlength_path = sys.argv[2]
    tax_cvrg_path = sys.argv[3]
    strain_reg_param = 1
    tax_fuzz_param = 1e-6
    seq_fuzz_param = 1e-6
    tax_uncertainty_param = 0.01

    info("Starting latent strain analysis.")
    info("seq_cvrg_path: %s", seq_cvrg_path)
    info("nlength_path: %s", nlength_path)
    info("tax_cvrg_path: %s", tax_cvrg_path)
    info("strain_reg_param: %f", strain_reg_param)
    info("tax_fuzz_param: %f", tax_fuzz_param)
    info("seq_fuzz_param: %f", seq_fuzz_param)
    info("tax_uncertainty_param: %f",tax_uncertainty_param)

    # Load data
    seq_cvrg  = pd.read_table(seq_cvrg_path,
                          names=['sample_id', 'sequence_id', 'tally'],
                          index_col=['sample_id', 'sequence_id'], squeeze=True
                         ).unstack('sequence_id', fill_value=0)
    seqlen = pd.read_table(nlength_path,
                            names=['contig_id', 'nlength'], index_col=['contig_id'], squeeze=True)
    seq_count = seq_cvrg.multiply(seqlen).round().astype(int)  # Scale to nucleotide counts
    tax_count = pd.read_table(tax_cvrg_path,
                          names=['sample_id', 'taxon_id', 'tally'],
                          index_col=['sample_id', 'taxon_id'], squeeze=True
                         ).unstack('taxon_id', fill_value=0)

    # Align tables
    seq_count = seq_count[seqlen.index]
    tax_count = tax_count.loc[seq_count.index]
    #
    # assert (seq_count.index == tax_count.index).all(), "Sequence and taxon table indices must be aligned."
    # assert (seqlen.index == seq_count.columns).all(), "Sequence and seqlen tables must be aligned."

    n_samples, g_seqs = seq_count.shape
    t_taxa = tax_count.shape[1]
    s_strains = 2 * t_taxa
    u_tax_counts = tax_count.sum(1).round().astype(int)
    v_seq_counts = seq_count.sum(1).round().astype(int)

    info("n_samples: %d", n_samples)
    info("g_seqs: %d", g_seqs)
    info("t_taxa: %d", t_taxa)
    info("s_strains: %d", s_strains)
    info("mean u_tax_counts: %f", u_tax_counts.mean())
    info("mean v_seq_counts: %f", v_seq_counts.mean())

    with pm.Model() as model:
        # Latent strain abundances
        alpha_raw = np.exp(-strain_reg_param * tt.arange(s_strains)/s_strains)
        alpha = alpha_raw / alpha_raw.sum()
        pi = pm.Dirichlet('pi', a=alpha, shape=(n_samples, s_strains))

        # Strains to taxa
        theta = pm.Dirichlet('theta',
                            a=np.ones(t_taxa) * tax_uncertainty_param / t_taxa,
                            shape=(s_strains, t_taxa))
        expect_tax_frac_unnorm = pi.dot(theta) + tax_fuzz_param
        obs_tax = pm.Multinomial('obs_tax',
                                n=u_tax_counts,
                                p=expect_tax_frac_unnorm,
                                shape=(n_samples, t_taxa),
                                observed=tax_count.values)

        # Strains to seqs
        phi = pm.Exponential('phi', lam=1, shape=(s_strains, g_seqs))
        expect_seq_frac_unnorm = pi.dot(phi) * seqlen + seq_fuzz_param


        obs_seq = pm.Multinomial('obs_seq',
                                n=v_seq_counts,
                                p=expect_seq_frac_unnorm,
                                shape=(n_samples, g_seqs),
                                observed=seq_count.values)

        # Save intermediate values:
        pm.Deterministic('expect_tax_frac', tt_simplex_normalize(expect_tax_frac_unnorm))
        pm.Deterministic('expect_seq_frac', tt_simplex_normalize(expect_seq_frac_unnorm))


    advi = pm.fit(int(1.5e5), model=model)
    info("Finished fitting model with ADVI.")
    advi_mean = advi.bij.rmap(advi.mean.eval())

    theta_est = model.theta.eval({model.theta_stickbreaking__: advi_mean['theta_stickbreaking__']})
    phi_est = model.phi.eval({model.phi_log__: advi_mean['phi_log__']})
    pi_est = model.pi.eval({model.pi_stickbreaking__: advi_mean['pi_stickbreaking__']})

    info("DONE")
