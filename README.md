# Strain Genomes from Metagenomes

## Purpose

This project is an exploration of a limited space of pipelines/models to
recover strain-resolved bacterial genomes from metagenomic reads.

The method will be based on a linear decomposition of "counts" into
components that represent individual strains and loadings that represent
the per-sample abundances of those strains.

## Notes about the general plan

### What is a strain?

TODO

### What is a count?

Sequence data always comes back to the non-negative integers (counts),
because you're inherently collecting discrete data.
Is it useful to exploit this information?
For one thing, this refines the "decomposition" framework, so I'd expect it to
be an improvement.

What discrete things might we be counting?

-   I think it might be a good idea to use kmer counts, which would forgo assembly,
    although this might result in too much data.
-   Another possibility would be assembly, gene calling, and read mapping to get
    counts, which would centralize
-   Yet another possibility is counting whole contigs
-   Or I could count bins (the least data intensive option)
-   _Or_, I could take some sort of hybrid approach where I decompose
    bins into "groups", and use just a subset of groups for downstream
    analyses (this is what I do in longev-mgen)

### How do I try this out?

I've got some data from longev-mgen that I could use.  This dataset also includes
matched 16S rRNA gene surveys, so I can use a cross-decomposition for potentially
more power.

### What do I do once I have latent components (strains)

One option would be to use these to filter a subset of reads and re-assemble.
Alternatively, I could just focus on getting gene presence absence.

## Notes about potentially useful literature

I'm seeing a whole bunch of references to statistical (ML) models that can be useful.
Key words include:

-   Sparsity
-   Topic modeling
-   Latent topics
-   Components
-   Counts
-   Latent Dirichlet allocation
-   SVN
-   Non-negative feature extraction

Primary among these:

-   [Latent semantic analysis](https://www.cs.cmu.edu/~jgc/publication/PublicationPDF/Sparse_Latent_Semantic_Analysis.pdf)
    (LSA) or topic modeling.  Basically treat the things I'm counting as words, the
    metagenomic libraries as documents,
    the genomes as latent topics
-   [Sparse Overcomplete Latent Variable Decomposition of Counts Data](https://paris.cs.illinois.edu/pubs/final_nips2007.pdf)
     Not sure how this is different from the other one

Some other stuff worth checking out:

-   https://scikit-learn.org/stable/modules/classes.html#module-sklearn.decomposition

## Notes on computationally efficient decompositions

It is important to note that the general problem I'm facing has
the potential to be very computationally intensive (both memory and CPU).
As such, it is worth thinking at large about what can make it more efficient.
Right now, I think that means thinking about (1) dimensional reduction
and (2) efficient data structures.
As key examples of (2), my data is almost always going to be sparse,
since not all taxa are at high abundance in all samples.

On this note, I'll be messing around with sparse matrices in scipy.
Right now, I think that entails loading my data in as a pandas DataFrame
((?) Is this too expensive?), and then transforming it into a DoK
sparse matrix representation, and subsequently a CSR or CSC format.

Using a sparse matrix is only really worthwhile, however, if the algorithm I'm using
can take advantage.
One example of such an algorithm is `sklearn.decomposition.TruncatedSVD`.
