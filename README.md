# Strain Genomes from Metagenomes

## Purpose

This project is an exploration of a limited space of pipelines/models to
recover strain-resolved bacterial genomes from metagenomic reads.

The method will be based on a linear decomposition of "counts" into
components that represent individual strains and loadings that represent
the per-sample abundances of those strains.

## What is a strain?

TODO

## What is a count?

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

## How do I try this out?

I've got some data from longev-mgen that I could use.  This dataset also includes
matched 16S rRNA gene surveys, so I can use a cross-decomposition for potentially
more power.

## What do I do once I have latent components (strains)

One option would be to use these to filter a subset of reads and re-assemble.
Alternatively, I could just focus on getting gene presence absence.
