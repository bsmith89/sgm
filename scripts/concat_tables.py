#!/usr/bin/env python3
"""Combine multiple tables with matching columns into one table.

Columns do not need to be in the same order, but they do need to have the same
header.

usage:

    concat_tables.py TABLE1 TABLE2 ... > OUTPUT

"""

import pandas as pd
import sys

def get_tables(paths):
    for p in paths:
        try:
            yield pd.read_table(p)
        except pd.io.common.EmptyDataError as err:
            print(p, file=sys.stderr)
            raise err

tables = get_tables(sys.argv[1:])
pd.concat(tables).to_csv(sys.stdout, sep='\t', index=False)
