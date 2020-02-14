import os
import json

import numpy as np
import networkx as nx
from Bio import SeqIO
import pysam

from .reads import get_covarying_sites
from .reads import obtain_superreads
from .graph import create_simple_reduced
from .graph import get_candidates
from .regression import perform_regression
from .regression import obtain_quasispecies


def write_json(path, data):
    if path:
        with open(path, 'w') as json_file:
            json.dump(data, json_file, indent=2)


def write_fasta(path, data):
    if path:
        SeqIO.write(data, path, 'fasta')


def write_csv(path, data):
    if path:
        data.to_csv(path)


def full_pipeline_io(
        input_bam, output_dir=None, output_fasta=None
        ):
    if output_dir and output_fasta:
        print(
            'WARNING: Ignoring fasta argument since output directory',
            'was specified'
        )
    alignment = pysam.AlignmentFile(input_bam, 'rb')
    covarying_sites, consensus = get_covarying_sites(alignment)
    superreads = obtain_superreads(alignment, covarying_sites)
    graph = create_simple_reduced(superreads)
    candidate_vacs, describing_superreads = get_candidates(graph)
    quasispecies_info = perform_regression(
        superreads, describing_superreads
    )
    all_quasispecies = obtain_quasispecies(
        quasispecies_info, superreads, consensus, covarying_sites
    )
    SeqIO.write(all_quasispecies, output_fasta, 'fasta')
