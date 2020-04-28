import json

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np
import pandas as pd
import pysam


characters = ['A', 'C', 'G', 'T', '-']


def single_read_count_data(read):
    sequence_length = np.array([
        cigar_tuple[1]
        for cigar_tuple in read.cigartuples
        if cigar_tuple[0] != 1 and cigar_tuple[0] != 4
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
    consensus_sequence = ''.join(nucleotide_counts['consensus'])
    consensus_record = SeqRecord(
        Seq(consensus_sequence),
        id='consensus',
        description=''
    )
    return covarying_sites[desired], consensus_record, nucleotide_counts


def read_reference_start_and_end(alignment, site_boundaries):
    read_information = pd.DataFrame(
        [
            (read.reference_start, read.reference_end)
            for read in alignment.fetch()
        ],
        columns=['reference_start', 'reference_end']
    )
    read_information['covarying_start'] = np.searchsorted(
        site_boundaries, read_information['reference_start']
    )
    read_information['covarying_end'] = np.searchsorted(
        site_boundaries, read_information['reference_end']
    )
    return read_information


def extract_label(query_name):
    return '-'.join(query_name.split('.')[:4]) if not '+' in query_name else 'AR'


def obtain_superreads(alignment, covarying_sites):
    read_information = read_reference_start_and_end(
        alignment, covarying_sites
    )
    read_groups = {}
    all_reads = list(alignment.fetch())
    for i, read in enumerate(all_reads):
        row = read_information.iloc[i, :]
        key = (row['covarying_start'], row['covarying_end'])
        if key in read_groups:
            read_groups[key].append(i)
        else:
            read_groups[key] = [i]
    all_superreads = []
    superread_index = 0
    for covarying_boundaries, read_group in read_groups.items():
        if covarying_boundaries[0] == covarying_boundaries[1]:
            continue
        superreads = {}
        for read_index in read_group:
            read = all_reads[read_index]
            label = extract_label(read.query_name)
            covarying_sites_in_read = covarying_sites[
                covarying_boundaries[0]: covarying_boundaries[1]
            ]
            value_at_covarying_sites = ''.join(
                [
                    read.query[triplet[0]].upper()
                    for triplet in read.get_aligned_pairs(True)
                    if triplet[1] in covarying_sites_in_read
                ]
            )
            has_ar = 1 if '+' in read.query_name else 0
            if value_at_covarying_sites in superreads:
                current_superread = superreads[value_at_covarying_sites]
                current_superread['weight'] += 1
                current_superread['ar'] += has_ar
                if not label in current_superread['composition']:
                     current_superread['composition'][label] = 0
                current_superread['composition'][label] += 1
                current_superread['read_names'].append(read.query_name)
            else:
                superreads[value_at_covarying_sites] = {
                    'weight': 1,
                    'ar': has_ar,
                    'composition': {label: 1},
                    'read_names': [read.query_name]
                }
        total_weight = sum([
            superread[1]['weight'] for superread in superreads.items()
        ])
        for vacs, info in superreads.items():
            all_superreads.append({
                'index': superread_index,
                'vacs': vacs,
                'weight': info['weight'],
                'frequency': info['weight']/total_weight,
                'ar': info['ar'],
                'ar_frequency': info['ar']/info['weight'],
                'cv_start': int(covarying_boundaries[0]),
                'cv_end': int(covarying_boundaries[1]),
                'composition': info['composition'],
                'read_names': info['read_names'],
                'discarded': False
            })
            superread_index += 1
    return all_superreads


def covarying_sites_io(
        bam_path, json_path, fasta_path, csv_path, threshold=.01
        ):
    alignment = pysam.AlignmentFile(bam_path, 'rb')
    covarying_sites, consensus, count_data = get_covarying_sites(
        alignment, threshold=threshold
    )
    covarying_sites_json = [int(site) for site in covarying_sites]
    with open(json_path, 'w') as json_file:
        json.dump(covarying_sites_json, json_file)
    SeqIO.write(consensus, fasta_path, 'fasta')
    count_data.to_csv(csv_path)


def superread_json_io(bam_path, covarying_path, superread_path):
    with open(covarying_path) as json_file:
        covarying_sites = np.array(json.load(json_file), dtype=np.int)
    alignment = pysam.AlignmentFile(bam_path, 'rb')

    superreads = obtain_superreads(alignment, covarying_sites)
    with open(superread_path, 'w') as json_file:
        json.dump(superreads, json_file, indent=2)


def superread_fasta_io(input_cvs, input_srdata, output_fasta,
        weight_filter=0, vacs_filter=0):
    with open(input_cvs) as json_file:
        cvs = json.load(json_file)
    with open(input_srdata) as json_file:
        srdata = json.load(json_file)
    outfile = open(output_fasta, 'w')
    n_cvs = len(cvs)
    for sr in srdata:
        heavy_enough = sr['weight'] > weight_filter
        long_enough = len(sr['vacs']) > vacs_filter
        should_write = heavy_enough and long_enough
        if should_write:
            outfile.write('>superread-%d_weight-%d\n' % (sr['index'], sr['weight']))
            seq = sr['cv_start']*'-' + sr['vacs'] + (n_cvs - sr['cv_end'])*'-'
            outfile.write(seq + '\n')
    outfile.close()
