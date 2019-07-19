import unittest
import os
import itertools as it

import pysam
import pandas as pd
import numpy as np
from Bio import SeqIO

from convex_qsr import ErrorCorrection
from convex_qsr import partial_covariation_test
from .mock import MockPysamAlignedSegment


data_dir = os.path.join('tests', 'data')


def get_numpy_fasta():
    fasta_path = os.path.join(data_dir, 'sorted.fasta')
    fasta = SeqIO.parse(fasta_path, 'fasta')
    return np.array([
        list(str(record.seq)) for record in fasta],
        dtype='<U1'
    )


def get_fasta_counts():
    numpy_fasta = get_numpy_fasta()
    counts = pd.DataFrame([
        pd.Series(numpy_fasta[:, i]).value_counts()
        for i in range(numpy_fasta.shape[1])
    ]).fillna(0)
    def zeros(character): return (counts[character] == 0).astype(np.int)
    counts['zero_cols'] = zeros('A') + zeros('C') + zeros('G') + zeros('T')
    counts['interesting'] = counts['zero_cols'] < 3
    return counts


def get_pairs_to_test(threshold):
    print('Generating matrix of pairs...')
    numpy_fasta = get_numpy_fasta()
    counts = get_fasta_counts()
    interesting_sites = counts.index[counts['interesting']]
    pairs = []
    all_pairs = it.combinations(interesting_sites, 2)
    for i, pair in enumerate(all_pairs):
        site_i, site_j = pair
        content_i = numpy_fasta[:, site_i] != '-'
        content_j = numpy_fasta[:, site_j] != '-'
        valid = content_i & content_j
        if valid.sum() > threshold:
            pairs.append((site_i, site_j))
    pairs_df = pd.DataFrame({
        'site_i': pd.Series([pair[0] for pair in pairs], dtype='int'),
        'site_j': pd.Series([pair[1] for pair in pairs], dtype='int')
    })
    return pairs_df


class TestErrorCorrection(unittest.TestCase):
    threshold = 20
    test_data = os.path.join(data_dir, 'sorted.bam')
    bam = pysam.AlignmentFile(test_data, 'rb')
    error_correction = None

    def setUp(self):
        if self.error_correction is None:
            self.__class__.error_correction = ErrorCorrection(self.bam)

    def test_single_read_count_data(self):
        segment = MockPysamAlignedSegment(
            'CTGATCGCTAACTA',  # read 8
            [(0, 4), (1, 1), (0, 6), (2, 1), (0, 3)],
            2
        )
        desired_sequence = np.array(list('CTGACGCTAA-CTA'), dtype='<U1')
        desired_positions = np.arange(2, 16)
        sequence, positions = ErrorCorrection.read_count_data(segment)

        sequences_agree = (sequence == desired_sequence).all()
        self.assertTrue(sequences_agree)

        positions_agree = (positions == desired_positions).all()
        self.assertTrue(positions_agree)

    def test_real_nucleotide_counts(self):
        bam_counts = self.error_correction.get_nucleotide_counts()
        desired_columns = ['A', 'C', 'G', 'T', 'interesting']

        fasta_counts = get_fasta_counts()
        fasta_subset = fasta_counts.loc[:, desired_columns]
        bam_subset = bam_counts.loc[:, desired_columns]
        fasta_equals_bam = fasta_subset == bam_subset
        self.assertTrue(fasta_equals_bam.all().all())

    def test_partial_covariation_test(self):
        pairs = self.error_correction.get_pairs()
        partial_covariation_test((self.bam.filename, pairs[:1000], 1))

    def test_full_covariation_test(self):
        self.error_correction.full_covariation_test()
        covarying_sites = self.error_correction.multiple_testing_correction()
        print(covarying_sites+1)

    def test_write_corrected_reads(self):
        corrected_bam_filename = 'corrected.bam'
        self.error_correction.write_corrected_reads(corrected_bam_filename)
        os.remove(corrected_bam_filename)
