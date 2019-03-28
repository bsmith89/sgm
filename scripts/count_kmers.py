#!/usr/bin/env python3

import sys
from collections import Counter
from functools import lru_cache
from Bio.SeqIO import parse
from tqdm import tqdm

COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def iter_kmers(string, k):
    for i in range(len(string) - k + 1):
        yield string[i:i+k]

def reverse_complement(string):
    return ''.join(COMPLEMENT[n] for n in reversed(string))

@lru_cache(maxsize=1024)
def regularize_kmer(kmer):
    '''Reverse complement kmer if superseded.'''
    rkmer = reverse_complement(kmer)
    if kmer <= rkmer:
        return kmer
    else:
        return rkmer

if __name__ == '__main__':
    for rec in tqdm(parse(sys.stdin, 'fasta')):
        kmer_counts = Counter(regularize_kmer(kmer)
                              for kmer
                              in iter_kmers(str(rec.seq), int(sys.argv[1])))
        for kmer in kmer_counts:
            print(rec.id, kmer, kmer_counts[kmer], sep='\t')
