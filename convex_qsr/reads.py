import json

import numpy as np
import pandas as pd
import pysam


characters = ['A', 'C', 'G', 'T', '-']


def single_read_count_data(read):
    sequence_length = np.array([
        cigar_tuple[1]
        for cigar_tuple in read.cigartuples
        if cigar_tuple[0] != 1
    ]).sum()
    first_position = read.reference_start
    last_position = first_position + sequence_length
    positions = np.arange(first_position, last_position)
    segments = []
    number_of_cigar_tuples = len(read.cigartuples)
    unaligned_sequence = read.query_alignment_sequence
    position = 0
    for i, cigar_tuple in enumerate(read.cigartuples):
        action = cigar_tuple[0]
        stride = cigar_tuple[1]
        match = action == 0
        insertion = action == 1
        deletion = action == 2
        if match:
            segments.append(
                unaligned_sequence[position: position + stride]
                )
            position += stride
        elif insertion:
            position += stride
        elif deletion:
            if len(segments) > 0 and i < number_of_cigar_tuples:
                segments.append(stride * '-')
    sequence = np.concatenate([list(segment) for segment in segments])
    return sequence, positions


def all_read_count_data(alignment):
    reference_length = alignment.header['SQ'][0]['LN']
    counts = np.zeros((reference_length, 5))
    for read in alignment.fetch():
        sequence, positions = single_read_count_data(read)
        for character_index, character in enumerate(characters):
            rows = positions[sequence == character]
            counts[rows, character_index] += 1
    return counts


def supplementary_info(row):
    result = row.loc[['A', 'C', 'G', 'T']] \
        .sort_values(ascending=False)
    result.index = ['c1', 'c2', 'c3', 'c4']
    return result.append(pd.Series(
        result.values/row['coverage'] if row['coverage'] else 0,
        index=['f1', 'f2', 'f3', 'f4']
    ))


def zeros(df, character):
    return (df[character] == 0).astype(np.int)


def site_table(alignment):
    counts = all_read_count_data(alignment)
    df = pd.DataFrame(counts, columns=characters)
    col_0s = zeros(df, 'A') + zeros(df, 'C') + zeros(df, 'G') + zeros(df, 'T')
    df['interesting'] = col_0s < 3
    df['nucleotide_max'] = df[['A', 'C', 'G', 'T']].max(axis=1)
    df['coverage'] = df[['A', 'C', 'G', 'T']].sum(axis=1)
    df['consensus'] = '-'
    for character in characters[:-1]:
        consensus_agreement = df['nucleotide_max'] == df[character]
        df.loc[consensus_agreement, 'consensus'] = character
    return pd.concat([df, df.apply(supplementary_info, axis=1)], axis=1)


def get_covarying_sites(alignment, threshold=.01, end_correction=10):
    nucleotide_counts = site_table(alignment)
    above_threshold = (
        nucleotide_counts
        .loc[:, ['f1', 'f2', 'f3', 'f4']] > threshold
    ).sum(axis=1)
    all_integers = np.arange(0, len(above_threshold))
    covarying_sites = all_integers[above_threshold > 1]
    number_of_sites = len(nucleotide_counts)
    after_head_correction = covarying_sites > end_correction
    final_site = number_of_sites - end_correction
    before_tail_correction = covarying_sites < final_site
    desired = after_head_correction & before_tail_correction
    return covarying_sites[desired]


def covarying_sites_io(bam_path, json_path):
    alignment = pysam.AlignmentFile(bam_path, 'rb')
    covarying_sites = get_covarying_sites(alignment)
    covarying_sites_json = [int(site) for site in covarying_sites]
    with open(json_path, 'w') as json_file:
        json.dump(covarying_sites_json, json_file)


