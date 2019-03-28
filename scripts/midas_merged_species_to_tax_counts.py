#!/usr/bin/env python3

import sys
import pandas as pd


if __name__ == '__main__':
    data_path = sys.argv[1]
    min_reads = int(sys.argv[2])

    data = pd.read_table(sys.argv[1], index_col='species_id')
    data = data[data.sum(1) > min_reads]
    data = data.rename(columns=lambda s: s[:-8])  # Drop '.m.midas' from merged table
    data.T.stack().to_csv(sys.stdout, sep='\t', header=False)
