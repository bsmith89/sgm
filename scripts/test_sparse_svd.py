#!/usr/bin/env python3

import sys

import pandas as pd
import numpy as np
import scipy as sp
from heapq import merge
from itertools import groupby
from sklearn.decomposition import TruncatedSVD

def parse_tallies(path, lib):
    with open(path) as handle:
        for line in handle:
            kmer, tally = line.split()
            tally = int(tally)
            yield kmer, lib, tally

def merge_tallies(*tallies):
    merge(*[sorted(t for t in tallies)])

def rows(merged_tallies):
    for kmer_group in groupby(merged_tallies, lambda x: x[0]):
        kmer = kmer_group[0]
        yield kmer, {lib: tally for _, lib, tally in kmer_group[1]}

if __name__ == '__main__':
    all_tallies = []
    libs = []
    for arg in sys.argv[1:]:
        lib, path = arg.strip().split(':')
        all_tallies.append(parse_tallies(path, lib))
        libs.append(lib)
    # Print header.
    print('', *libs, sep='\t')
    for kmer, by_lib in rows(merge(*all_tallies)):
        # Print index (kmer).
        print(kmer, end='')
        for l in libs:
            # Print elements (tallies).
            if l in by_lib:
                print('\t{}'.format(by_lib[l]), end='')
            else:
                print('\t0', end='')
        # Finish the line.
        print()

