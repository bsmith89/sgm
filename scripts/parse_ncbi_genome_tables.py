#!/usr/bin/env python
# coding: utf-8

import pandas as pd

rename_map = { 'Organism Name': 'taxon_name'
                    , 'Organism Groups': 'taxonomy_string'
                    , 'Strain': 'strain_name'
#                    , 'Level': 'assembly_completeness'
#                    , 'Size(Mb)': 'total_length'
#                    , 'GC%': 'gc_percent'
                    , 'Replicons': 'replicons_string'
#                    , 'Scaffolds': 'scaffold_tally'
#                    , 'CDS': 'cds_tally'
#                    , 'Release Date': 'release_date'
#                    , 'GenBank FTP': 'genbank_ftp_url'
                    , 'RefSeq FTP': 'refseq_ftp_url'
#                    , 'Genes': 'gene_tally
#                    , 'Host': 'host_name'
#                    , 'Modify Date': 'modify_date'
                    }


def taxon_name_to_id(name):
    'Transform name from NCBI into a normalized genome_id'
    out = name
    words = out.split(' ')
    # Filter out non-alpha characters from genus word.
    words[0] = ''.join(filter(str.isalpha, words[0]))
    # Replace genus with first letter.
    words[0] = words[0][0]
    # Drop the part after a '=' word (synonym?)
    if '=' in words:
        words = words[:words.index('=')]
    # Remove uneccessary modifier words.
    words = filter(lambda s: s not in ['str.', 'substr.', 'subsp.', 'bv.', 'biovar'], words)
    out = '_'.join(words)
    # Replace with better characters
    for char in [':', '-', '/', '=', '.']:
        out = out.replace(char, '_')
    for char in ['(', ')', "'"]:
        out = out.replace(char, '')

    # Check only legal characters
    for char in out:
        assert char.isalnum() or char == '_', (name, out)

    return out


def parse_replicon(r):
    r = r.strip()
    name, accessions = r.split(':')
    if name.startswith('chromosome'):
        replicon_type = 'chromosome'
    elif name.startswith('plasmid'):
        replicon_type = 'plasmid'
    if '/' in accessions:
        refseq_id, genbank_id = accessions.split('/')
    else:
        genbank_id = accessions
        refseq_id = ''
    return name, replicon_type, refseq_id, genbank_id


def flatten_replicons(x):
    # replicons string is "<chromosome name>" or "plasmid" followed by ":<refseq>/<genbank>"
    genome_id = x.genome_id
    replicon = [parse_replicon(r) for r in x.replicons_string.split(';')]
    for r in replicon:
        yield tuple([genome_id]) + r

if __name__ == "__main__":
    data = pd.read_csv('meta/ncbi_genomes.csv').rename(columns={'#Organism Name': 'Organism Name'})
    data = data.rename(columns=rename_map)

    data['genome_id'] = data.taxon_name.map(taxon_name_to_id)
    assert data.genome_id.is_unique

    replicon = []
    for _, g in data.iterrows():
        for r in flatten_replicons(g):
            replicon.append(r)

    replicon = pd.DataFrame(replicon, columns=['genome_id', 'replicon_name', 'replicon_type', 'refseq_id', 'genbank_id'])

    # Deal with some special cleanup cases.
    # Two plasmids have the same name.
    replicon.loc[replicon.genbank_id == 'CP018754.1', 'replicon_name'] += '_1'
    replicon.loc[replicon.genbank_id == 'CP018755.1', 'replicon_name'] += '_2'
    # Two chromoses have the same genbank_id, but one, despite having a different refseq_id, doesn't really exist.
    replicon = replicon.drop(replicon[replicon.refseq_id == 'NZ_CP015586.1'].index)

    # Produce unique replicon_ids.
    replicon['replicon_id'] = replicon.groupby('genome_id').genome_id.apply(lambda x: x + ['_' + str(i) for i in range(len(x))])

    # Check no duplicates for three different unique 'columns'.
    for ser in [replicon.replicon_id,
                (replicon.genome_id + replicon.replicon_name),
                replicon.refseq_id,
                replicon.genbank_id]:
        assert replicon[ser.duplicated(False)].empty

    genome = data[['genome_id', 'taxon_name', 'taxonomy_string', 'strain_name', 'refseq_ftp_url']]
    for col in ['genome_id', 'taxon_name', 'refseq_ftp_url']:
        assert genome[col].is_unique

    # Output tables.
    genome[['genome_id', 'taxon_name',
            'taxonomy_string', 'strain_name',
            'refseq_ftp_url']].to_csv('data/genome.tsv', sep='\t', index=False)
    replicon[['replicon_id', 'genome_id',
            'replicon_name', 'replicon_type',
            'refseq_id', 'genbank_id']].to_csv('data/replicon.tsv', sep='\t', index=False)
