#!/usr/bin/env python3

import sys


def chunk_depth(handle):
    """Split depth file lines into tables by the contig_id."""
    chunk = []
    for line in handle:
        contig_id, position, depth = line.split()
        if not chunk:
            chunk_contig_id = contig_id
        line_contig_id = contig_id
        if line_contig_id == chunk_contig_id:
            chunk.append((contig_id, int(position), int(depth)))
        else:
            yield chunk
            chunk = [(contig_id, int(position), int(depth))]
            chunk_contig_id = line_contig_id
    yield chunk


def main():
    depth_path, length_path, float_fmt = sys.argv[1:]
    print('contig_id', 'coverage', sep='\t', file=sys.stdout)

    with open(depth_path) as depth_handle, \
            open(length_path) as length_handle:
        for depth_table in chunk_depth(depth_handle):
            contig_id = depth_table[0][0]
            contig_id_b = None
            # Scan for the current contig in the metadata
            # (This depends on contigs being in the same order in both files.)
            while contig_id_b != contig_id:
                contig_id_b, length = next(length_handle).split()
            length = int(length)
            # Transpose the table, and extract a list of positions and depths.
            _, positions, depths = zip(*depth_table)
            coverage = sum(depths) / length
            # num_positions = len(positions)
            # hit_fraction = num_positions / length
            # left = positions[0] - 1
            # right = positions[-1]
            # span_length = (right - left)
            # contiguous = int(span_length == num_positions)
            print(contig_id, float_fmt % coverage, sep='\t', file=sys.stdout)


if __name__ == '__main__':
    main()
