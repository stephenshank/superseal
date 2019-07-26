import json

import numpy as np
import networkx as nx
from Bio import SeqIO
import pysam

from .error_correction import ErrorCorrection
from .mapped_reads import MappedReads
from .read_graph import SuperReadGraph
from .regression import perform_regression


def error_correction_io(
        input_bam, output_bam, output_json=None, output_fasta=None,
        end_correction=None
        ):
    alignment = pysam.AlignmentFile(input_bam, 'rb')
    error_correction = ErrorCorrection(
        alignment, end_correction=end_correction
    )
    error_correction.write_corrected_reads(output_bam)
    pysam.index(output_bam)

    if output_json:
        with open(output_json, 'w') as json_file:
            json.dump(
                [int(i) for i in error_correction.covarying_sites],
                json_file
            )

    if output_fasta:
        consensus_record = error_correction.consensus()
        SeqIO.write(consensus_record, output_fasta, 'fasta')


def read_graph_io(
        input_bam, input_json, input_consensus,
        output_full, output_cvs, output_describing,
        output_graph, output_candidates
        ):
    with open(input_json) as json_file:
        covarying_sites = np.array(json.load(json_file), dtype=np.int)
    consensus = SeqIO.read(input_consensus, 'fasta')
    corrected_reads = MappedReads(input_bam, 'rb')

    superread_graph = SuperReadGraph(corrected_reads, covarying_sites)
    superread_graph.obtain_superreads()
    superread_graph.create_full()

    full_superreads = superread_graph.get_full_superreads(consensus)
    SeqIO.write(full_superreads, output_full, 'fasta')

    superread_graph.reduce()
    cvs_superreads = superread_graph.get_cvs_superreads()
    SeqIO.write(cvs_superreads, output_cvs, 'fasta')

    superread_json = nx.node_link_data(superread_graph.superread_graph)
    with open(output_graph, 'w') as json_file:
        json.dump(superread_json, json_file, indent=2)

    superread_records, describing_superreads = \
        superread_graph.get_candidate_quasispecies(consensus)
    with open(output_describing, 'w') as json_file:
        json.dump(describing_superreads, json_file, indent=2)
    SeqIO.write(superread_records, output_candidates, 'fasta')


def regression_io(
        input_superreads, input_describing_json,
        input_candidates_fasta, output_fasta
        ):
    with open(input_superreads) as superread_file:
        superreads = json.load(superread_file)['nodes']

    with open(input_describing_json) as candidates_file:
        describing = json.load(candidates_file)

    quasispecies_info = perform_regression(superreads, describing)
    candidates = list(SeqIO.parse(input_candidates_fasta, 'fasta'))
    all_quasispecies = []
    for i, quasispecies in enumerate(quasispecies_info):
        record = candidates[quasispecies['index']]
        record_info = (i + 1, quasispecies['frequency'])
        record.id = 'quasispecies-%d_frequency-%.5f' % record_info
        record.description = ''
        all_quasispecies.append(record)

    SeqIO.write(all_quasispecies, output_fasta, 'fasta')
