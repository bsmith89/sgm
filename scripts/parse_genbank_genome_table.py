#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import sys

def construct_genome_id(x):
    words = x.organism_name.split(' ')

    # Drop useless words
    words = [w for w in words if w not in ['sp.']]
    genus, species, *_ = words

    # Collect strain name
    strain = x.infraspecific_name[len('strain='):]

    out = f'{genus}_{species}_{strain}'

    # Replace spacing chars
    for char in [' ', '/', '-', '.', '=']:
        out = out.replace(char, '_')

    # Replace some special cases:
    out = out.replace('+', 'pos')
    out = out.replace('*', 'star')

    # Drop other chars
    for char in ['(', ')', '[', ']', '#', ':', "'", ';', ',']:
        out = out.replace(char, '')

    while '__' in out:
        out = out.replace('__', '_')

    # Check chars
    for char in out:
        assert char.isalnum() or char == '_', f'Non-alphanumeric char found in {x.organism_name} {x.infraspecific_name} ({out})'
    return out


if __name__ == "__main__":
    data = (pd.read_table(sys.argv[1], skiprows=1)
            .rename(columns={'# assembly_accession': 'assembly_accession'})
            [lambda x: (x.assembly_level.isin(['Complete Genome'])) &
                        (x.genome_rep.isin(['Full'])) &
                        (x.infraspecific_name.str.startswith('strain=') &
                        (x.excluded_from_refseq.isna()))
            ]
        )

    data['genome_id'] = data.apply(construct_genome_id, axis=1)
    data['ftp_basename_stem'] = data.ftp_path.apply(lambda x: x.rsplit('/', 1)[-1])
    data['ftp_stem'] = data.ftp_path + '/' + data.ftp_basename_stem
    data['seq_rel_date'] = pd.to_datetime(data.seq_rel_date)
    data = data.sort_values('seq_rel_date', ascending=False).drop_duplicates(subset=['genome_id'])
    assert data.genome_id.is_unique

    data[['genome_id', 'organism_name',
        'infraspecific_name',
        'ftp_stem']].to_csv(sys.stdout, sep='\t', index=False)
