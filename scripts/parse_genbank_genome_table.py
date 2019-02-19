#!/usr/bin/env python3
# coding: utf-8

import pandas as pd
import sys
from tqdm import trange
import os

def construct_genome_id(x):
    # Drop useless words
    useless = ['sp.', 'subsp.', 'str.',
               'substr.', 'pv.', 'cf.',
               'serovar', 'biovar']
    species_words = [w
                     for w in x.organism_name.split(' ')
                     if w not in useless]
    genus, species, *_ = species_words

    # Collect strain name
    strain = x.infraspecific_name[len('strain='):]
    strain_words = [w
                    for w in strain.split(' ')
                    if w not in useless]
    strain = '_'.join(strain_words)

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

def fetch_replicon_data(x, dirname):
    path = os.path.join(dirname, x['ftp_basename_stem'] + '_assembly_report.txt')
    with open(path) as handle:
        for line in handle:
            if line.startswith('# Sequence-Name\tSequence-Role'):
                break
        d = pd.read_table(handle, names=['Sequence-Name', 'Sequence-Role',
                                       'Assigned-Molecule', 'Assigned-Molecule-Location/Type',
                                       'GenBank-Accn', 'Relationship', 'RefSeq-Accn',
                                       'Assembly-Unit', 'Sequence-Length', 'UCSC-style-name'],
                          dtype=str)
    d = d.rename(columns={'Assigned-Molecule-Location/Type': 'replicon_type',
                          'GenBank-Accn': 'genbank_id',
                          'Sequence-Length': 'size',
                          'Sequence-Name': 'replicon_name'})
    d.replicon_type = d.replicon_type.fillna('unknown')
    d['size'] = d['size'].astype(int)
    d.replicon_type = d.replicon_type.str.lower()
    # Add column with a incremental value for each replicon type
    d['replicon_index_by_type'] = (d.sort_values('genbank_id')  # Standardize order
                                    .groupby('replicon_type')  # Different incrementer per type
                                    .apply(lambda x: pd.DataFrame({'iter': range(len(x))},
                                                                  index=x.index)).iter + 1
                                  )
    d['replicon_anon_name'] = d['replicon_type'] + '_' + d['replicon_index_by_type'].astype('str')
    # TODO: Improve replicon names (currently a single chromosome will be named
    # "chromosome_1" despite there not being a chromosome_2.
    d['replicon_name'] = d.replicon_name.where(((d.replicon_name != 'ANONYMOUS') &
                                                (d.replicon_name.notna()) &
                                                (~d.replicon_name.str.startswith('unnamed'))),
                                               d.replicon_anon_name)
    d['replicon_id'] = x['genome_id'] + '_' + d['replicon_name']
    d['genome_id'] = x.genome_id
    return d[['replicon_id', 'genome_id', 'genbank_id', 'replicon_name', 'replicon_type', 'size']]


if __name__ == "__main__":
    print('Parsing genome table.', file=sys.stderr)
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

    print(f'Parsing assembly reports.', file=sys.stderr)
    replicon = []
    for i in trange(len(data)):
        rec = data.iloc[i]
        replicon_table = fetch_replicon_data(rec, sys.argv[2])
        replicon.append(replicon_table)
    replicon = pd.concat(replicon)

    data[['genome_id', 'organism_name',
        'infraspecific_name',
        'ftp_stem']].to_csv(sys.argv[3], sep='\t', index=False)
    replicon[['replicon_id', 'genome_id', 'genbank_id',
              'replicon_type', 'size']].to_csv(sys.argv[4], sep='\t', index=False)
