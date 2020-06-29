"Module for dealing with aligned reads."
import json
from collections import Counter

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import numpy as np
import pandas as pd
import pysam
import vcf

from .io import write_json


characters = ['A', 'C', 'G', 'T', '-']


def single_read_count_data(read):
    "Handle a single read."
    segments = []
    number_of_cigar_tuples = len(read.cigartuples)
    unaligned_sequence = read.query_alignment_sequence
    read_position = 0
    reference_position = read.reference_start
    positions = []
    for i, cigar_tuple in enumerate(read.cigartuples):
        action = cigar_tuple[0]
        stride = cigar_tuple[1]
        match = action == 0
        insertion = action == 1
        deletion = action == 2
        softclip = action == 4
        hardclip = action == 5

        if match:
            segments.append(
                unaligned_sequence[read_position: read_position + stride]
            )
            read_position += stride
            positions.append(
                np.arange(reference_position, reference_position + stride)
            )
            reference_position += stride
        elif insertion:
            read_position += stride
        elif deletion:
            if len(segments) > 0 and i < number_of_cigar_tuples:
                segments.append(stride * '-')
                positions.append(
                    np.arange(reference_position, reference_position + stride)
                )
                reference_position += stride
        elif softclip:
            read_position += stride
        elif hardclip:
            pass
    sequence = np.concatenate([list(segment) for segment in segments])
    positions = np.concatenate(positions)
    return sequence, positions


def all_read_count_data(alignment):
    "Process all reads."
    reference_length = alignment.header['SQ'][0]['LN']
    counts = np.zeros((reference_length, 5))
    for read in alignment.fetch():
        sequence, positions = single_read_count_data(read)
        for character_index, character in enumerate(characters):
            rows = positions[sequence == character]
            counts[rows, character_index] += 1
    return counts


def supplementary_info(row):
    "Get additional information derived from counts."
    result = row.loc[['A', 'C', 'G', 'T']] \
        .sort_values(ascending=False)
    result.index = ['c1', 'c2', 'c3', 'c4']
    return result.append(pd.Series(
        result.values/row['coverage'] if row['coverage'] else 0,
        index=['f1', 'f2', 'f3', 'f4']
    ))


def zeros(site_data, character):
    "Count zeros associated with a nucleotide."
    return (site_data[character] == 0).astype(np.int)


def site_table(alignment):
    "Get extra information associated with each site."
    counts = all_read_count_data(alignment)
    count_data = pd.DataFrame(counts, columns=characters)
    col_0s = zeros(count_data, 'A') + zeros(count_data, 'C') + \
        zeros(count_data, 'G') + zeros(count_data, 'T')
    count_data['interesting'] = col_0s < 3
    count_data['nucleotide_max'] = count_data[['A', 'C', 'G', 'T']].max(axis=1)
    count_data['coverage'] = count_data[['A', 'C', 'G', 'T']].sum(axis=1)
    count_data['consensus'] = '-'
    for character in characters[:-1]:
        consensus_agreement = count_data['nucleotide_max'] == count_data[character]
        count_data.loc[consensus_agreement, 'consensus'] = character
    return pd.concat([count_data, count_data.apply(supplementary_info, axis=1)], axis=1)


def get_covarying_sites(
        alignment, threshold=.01, end_correction=10,
        minimum_coverage=100):
    "Determine covarying sites from alignment."
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
    coverage = nucleotide_counts.loc[covarying_sites, 'coverage']
    well_covered = coverage > minimum_coverage
    desired = after_head_correction & before_tail_correction & well_covered
    consensus_sequence = ''.join(nucleotide_counts['consensus'])
    consensus_record = SeqRecord(
        Seq(consensus_sequence),
        id='consensus',
        description=''
    )
    return covarying_sites[desired], consensus_record, nucleotide_counts


def read_reference_start_and_end(alignment, site_boundaries):
    "Locate what covarying sites are covered by each read."
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
    "Label reads with known truth."
    if not '+' in query_name:
        return '-'.join(query_name.split('.')[:4])
    return 'AR'


def obtain_superreads(alignment, covarying_sites):
    "Get superreads from reads."
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
                    read.seq[triplet[0]].upper()
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
    "IO function for covarying sites."
    alignment = pysam.AlignmentFile(bam_path, 'rb')
    covarying_sites, consensus, count_data = get_covarying_sites(
        alignment, threshold=threshold
    )
    covarying_sites_json = [int(site) for site in covarying_sites]
    if json_path:
        with open(json_path, 'w') as json_file:
            json.dump(covarying_sites_json, json_file)
    if fasta_path:
        SeqIO.write(consensus, fasta_path, 'fasta')
    if csv_path:
        count_data.to_csv(csv_path)


def covariation_input(covarying_path):
    "Helper function for parsing covariation input."
    if covarying_path.split('.')[-1] == 'json':
        with open(covarying_path) as json_file:
            covarying_sites = np.array(json.load(json_file), dtype=np.int)
        return covarying_sites
    vcf_reader = vcf.Reader(filename=covarying_path)
    position_counter = Counter()
    for variant in vcf_reader:
        position_counter[variant.POS - 1] += 1
    unique_ordered_sites = sorted([
        site
        for site, count in position_counter
        if count > 1
    ])
    return np.array(unique_ordered_sites)


def superread_json_io(bam_path, covarying_path, superread_path):
    "IO function for creating superreads."
    covarying_sites = covariation_input(covarying_path)
    alignment = pysam.AlignmentFile(bam_path, 'rb')

    superreads = obtain_superreads(alignment, covarying_sites)
    with open(superread_path, 'w') as json_file:
        json.dump(superreads, json_file, indent=2)


def superread_fasta_io(
        input_cvs, input_superreaddata, output_fasta,
        weight_filter=0, vacs_filter=0):
    "IO function for creating superread FASTA file."
    cvs = covariation_input(input_cvs)
    with open(input_cvs) as json_file:
        cvs = json.load(json_file)
    with open(input_superreaddata) as json_file:
        superreaddata = json.load(json_file)
    outfile = open(output_fasta, 'w')
    n_cvs = len(cvs)
    for superread in superreaddata:
        heavy_enough = superread['weight'] > weight_filter
        long_enough = len(superread['vacs']) > vacs_filter
        should_write = heavy_enough and long_enough
        if should_write:
            outfile.write(
                '>superread-%d_weight-%d\n' %
                (superread['index'], superread['weight'])
            )
            seq = superread['cv_start']*'-' + \
                superread['vacs'] + \
                (n_cvs - superread['cv_end'])*'-'
            outfile.write(seq + '\n')
    outfile.close()


def resolvable_regions(superreads):
    "Determine resolvable regions from superreads."
    number_of_covarying_sites = max([
        sr['cv_end'] for sr in superreads
    ])
    pair_counts = np.zeros(number_of_covarying_sites - 1)
    for superread in superreads:
        pair_counts[
            superread['cv_start']: superread['cv_end'] - 1
        ] += superread['weight']
    insufficient_coverage = list(np.where(pair_counts == 0)[0])
    starts = [0] + insufficient_coverage
    stops = insufficient_coverage + [number_of_covarying_sites]
    regions = [
        {'start': int(start), 'stop': int(stop)}
        for start, stop in zip(starts, stops)
    ]
    return {
        'pair_counts': [int(pair) for pair in pair_counts],
        'regions': regions
    }


def resolvable_regions_io(input_srdata, output_json):
    "IO function for determining resolvable regions."
    with open(input_srdata) as json_file:
        superreads = json.load(json_file)
    regions = resolvable_regions(superreads)
    write_json(output_json, regions)
