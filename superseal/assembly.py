"Assembling quasispecies from superreads."
import json

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from .reads import characters
from .io import read_json
from .io import write_json


character_to_index_map = {char: i for i, char in enumerate(characters)}
index_to_character_map = characters


class Scaffold:
    "Scaffolds of intrahost variants."

    def __init__(self, superreads, members=None):
        self.members = members or {}
        self.superreads = superreads
        self.number_of_covarying_sites = max(sr['cv_end'] for sr in superreads)
        self.coverage = np.zeros(self.number_of_covarying_sites, dtype=np.int)

    def merge_node(self, superread):
        "Add node to this scaffold."
        if superread['filtered_index'] in self.members:
            return
        indices = np.arange(superread['cv_start'], superread['cv_end'])
        self.coverage[indices] += superread['weight']
        self.members[superread['filtered_index']] = True

    def merge_edge(self, superread_i, superread_j):
        "Add an edge to this scaffold."
        self.merge_node(superread_i)
        self.merge_node(superread_j)

    def merge_scaffold(self, scaffold):
        "Add another scaffold to this scaffold."
        for i in scaffold.members.keys():
            self.members[i] = True
        self.coverage += scaffold.coverage

    def check_membership(self, superread):
        "Check if a superread is a member of this scaffold."
        return superread['filtered_index'] in self.members

    def check_membership_by_index(self, superread_index):
        "Check if a superread is a member of this scaffold by its index."
        return superread_index in self.members

    def is_covered(self):
        "Check if this scaffold fully covers the region."
        return np.all(self.coverage > 0)

    def member_superreads(self):
        "Get all members of this scaffold."
        return [self.superreads[i] for i in self.members.keys()]

    def leftmost(self):
        "Find leftmost member of this scaffold."
        return min(self.member_superreads(), key=lambda sr: sr['cv_start'])

    def rightmost(self):
        "Find rightmost member of this scaffold."
        return max(self.member_superreads(), key=lambda sr: sr['cv_end'])

    def counts(self):
        "Get coverage counts."
        counts = np.zeros((self.number_of_covarying_sites, 5))
        for superread_index in self.members.keys():
            superread = self.superreads[superread_index]
            weight = superread['weight']
            for vac_index, vac in enumerate(superread['vacs']):
                character_index = character_to_index_map[vac]
                site_index = superread['cv_start'] + vac_index
                counts[site_index, character_index] += weight
        return counts

    def consensus(self):
        "Get consensus sequence of this scaffold."
        counts = self.counts()
        consensus = np.array(
            self.number_of_covarying_sites * ['N'], dtype='<U1'
        )
        for site_index in np.arange(self.number_of_covarying_sites):
            character_index = np.argmax(counts[site_index, :])
            character = index_to_character_map[character_index]
            consensus[site_index] = character
        return consensus

    def extremities(self):
        "Get extremal members of this scaffold."
        left = self.leftmost()
        right = self.rightmost()
        left_initiates = left['cv_start'] == 0
        right_terminates = right['cv_end'] == self.number_of_covarying_sites
        if left_initiates and not right_terminates:
            return {'left': None, 'right': right}
        elif right_terminates and not left_initiates:
            return {'left': left, 'right': None}
        elif left_initiates and right_terminates:
            return {'left': None, 'right': None}
        else:
            return {'left': left, 'right': right}


def superread_filter(minimum_weight, region):
    "Get a filter function for superreads based on weight and region."
    def the_actual_filter(superread):
        heavy_enough = superread['weight'] >= minimum_weight
        starts_after = superread['cv_start'] >= region['start']
        ends_before = superread['cv_end'] <= region['stop']
        desired = heavy_enough and starts_after and ends_before
        return desired
    return the_actual_filter


def filter_superreads(superreads, minimum_weight, region):
    "Filter superreads based on weight and region."
    filtered_superreads = list(filter(
        superread_filter(minimum_weight, region),
        superreads
    ))
    number_of_covarying_sites = 0
    for i, superread in enumerate(filtered_superreads):
        superread['filtered_index'] = i
        number_of_covarying_sites = max(
            number_of_covarying_sites, superread['cv_end']
        )
    return filtered_superreads, number_of_covarying_sites


def check_compatability(superread_i, superread_j, minimum_overlap=2):
    "Check whether two superreads are compatible."
    if not 'index' in superread_i or not 'index' in superread_j:
        return (False, 0)
    if superread_i['index'] == superread_j['index']:
        return (False, 0)
    i_cv_start = superread_i['cv_start']
    i_cv_end = superread_i['cv_end']
    j_cv_start = superread_j['cv_start']
    j_cv_end = superread_j['cv_end']
    start_before_start = i_cv_start < j_cv_start
    start_before_end = j_cv_start < i_cv_end
    end_before_end = i_cv_end < j_cv_end
    if start_before_start and start_before_end and end_before_end:
        cv_start = max(i_cv_start, j_cv_start)
        cv_end = min(i_cv_end, j_cv_end)
        delta = cv_end - cv_start
        i_start = cv_start - i_cv_start
        i_end = i_start + delta
        j_start = cv_start - j_cv_start
        j_end = j_start + delta
        i_sequence = superread_i['vacs'][i_start: i_end]
        j_sequence = superread_j['vacs'][j_start: j_end]
        agree_on_overlap = i_sequence == j_sequence
        overlap = len(i_sequence)
        weight = min(superread_i['weight'], superread_j['weight'])
        support = overlap * weight
        long_enough = overlap >= minimum_overlap
        compatible = agree_on_overlap and long_enough
        edge_data = {
            'weight': weight,
            'overlap': overlap,
            'support': support
        }
        return (compatible, edge_data)
    return (False, {})


def get_edge_list(superreads):
    "Get edge list from superreads."
    edge_list = []
    for i, superread_i in enumerate(superreads):
        for j, superread_j in enumerate(superreads):
            should_include_edge, edge_data = check_compatability(
                superread_i, superread_j
            )
            if should_include_edge:
                edge_list.append({
                    'i': i,
                    'j': j,
                    'available': True,
                    **edge_data
                })
    return pd.DataFrame(edge_list).sort_values(by='overlap', ascending=False)


def check_integrity(scaffolds):
    "Ensure scaffolds are disjoint."
    for i, scaffold_i in enumerate(scaffolds):
        for member in scaffold_i.members.keys():
            for scaffold_j in scaffolds[i+1:]:
                assert not scaffold_j.check_membership_by_index(member)


def scaffold_descriptions(scaffolds, all_superreads):
    "Get superread information about the makeup of all scaffolds."
    superreads = scaffolds[0].superreads
    number_of_covarying_sites = scaffolds[0].number_of_covarying_sites
    utilized = set([
        member
        for scaffold in scaffolds
        for member in list(scaffold.members.keys())
    ])
    unused = [i for i in range(len(all_superreads)) if not i in utilized]
    return {
        'describing_superreads': [
            sorted([int(member) for member in scaffold.members.keys()])
            for scaffold in scaffolds
        ],
        'original_indices': [
            sorted([
                int(superreads[member]['index'])
                for member in scaffold.members.keys()
            ])
            for scaffold in scaffolds
        ],
        'coverage': [
            [int(i) for i in scaffold.coverage] for scaffold in scaffolds
        ],
        'unused': unused,
        'full_coverage': [
            bool(np.all(scaffold.coverage)) for scaffold in scaffolds
        ],
        'consensus': [''.join(scaffold.consensus()) for scaffold in scaffolds],
        'number_of_covarying_sites': number_of_covarying_sites
    }


def invalidate_edges(edges, sr_index):
    "Mark edges that may no longer be used."
    uses_node = (edges['i'] == sr_index) | (edges['j'] == sr_index)
    edges.loc[uses_node, 'available'] = False


def get_initial_node(superreads, node_hash):
    "Get an initial seed node."
    superread = max(
        [sr for sr in superreads if node_hash[sr['filtered_index']]],
        key=lambda sr: sr['weight']
    )
    node_hash[superread['filtered_index']] = False
    return superread


def get_nodes_for_extension(edges, scaffold, metric='support'):
    "Get nodes that can extended."
    extremal_nodes = scaffold.extremities()
    available_edges = edges.loc[edges.available, :]
    nodes = []
    if not extremal_nodes['left'] is None:
        left_index = extremal_nodes['left']['filtered_index']
        extends_left = available_edges['j'] == left_index
        left_extending_edges = available_edges.loc[extends_left, metric]
        if len(left_extending_edges) > 0:
            best_edge = left_extending_edges.nlargest(1).index[0]
            nodes.append(available_edges.loc[best_edge, 'i'])
    if not extremal_nodes['right'] is None:
        right_index = extremal_nodes['right']['filtered_index']
        extends_right = available_edges['i'] == right_index
        right_extending_edges = available_edges.loc[extends_right, metric]
        if len(right_extending_edges) > 0:
            best_edge = right_extending_edges.nlargest(1).index[0]
            nodes.append(available_edges.loc[best_edge, 'j'])
    return nodes


def assemble_single_region(
        all_superreads, region, minimum_weight=5, max_qs=2,
        verbose=True):
    "Assemble an individual region."
    superreads, number_of_covarying_sites = filter_superreads(
        all_superreads, minimum_weight, region
    )
    edges = get_edge_list(superreads)
    node_hash = {i: True for i in range(len(superreads))}
    scaffolds = []
    stop_early = False
    for i in range(max_qs):
        if verbose:
            print('Scaffold %d of %d...\n' % (i+1, max_qs), end='')
        scaffold = Scaffold(superreads)
        initial_superread = get_initial_node(superreads, node_hash)
        if verbose:
            message = 'Seed superread for quasispecies %d: %d'
            data = (i+1, initial_superread['index'])
            print(message % data)
            print('Percent of covarying sites covered:', end='')
        scaffold.merge_node(initial_superread)

        # Expansion
        while not scaffold.is_covered():
            new_nodes = get_nodes_for_extension(edges, scaffold)
            for node_index in new_nodes:
                scaffold.merge_node(superreads[node_index])
                node_hash[node_index] = False
            if verbose:
                message = '%.2f '
                data = np.sum(scaffold.coverage > 0) / \
                    number_of_covarying_sites
                print(message % data, end='')

            if len(new_nodes) == 0 and not scaffold.is_covered():
                stop_early = True
                break

        if stop_early:
            if verbose:
                print('Stopping early... unable to fully cover!')
            return scaffold_descriptions(scaffolds[:-1], False)

        # Absorption
        consensus = scaffold.consensus()
        number_absorbed = 0
        for superread in [sr for sr in superreads if node_hash[sr['filtered_index']]]:
            relevant_consensus = ''.join(consensus[
                superread['cv_start']: superread['cv_end']
            ])
            if superread['vacs'] == relevant_consensus:
                number_absorbed += 1
                scaffold.merge_node(superread)
                node_hash[superread['filtered_index']] = False
        if verbose:
            remaining = np.sum(np.array(list(node_hash.values())))
            print('...done! Absorbed %d superreads.' % number_absorbed)
            print(remaining, 'remaining.')
        scaffolds.append(scaffold)
    print('Obtained superreads that describe %d quasispecies!!' % max_qs)
    return scaffold_descriptions(scaffolds, all_superreads)


def assemble(all_superreads, resolution, **kwargs):
    "Assemble all regions."
    all_descriptions = []
    for region in resolution['regions']:
        all_descriptions.append(
            assemble_single_region(all_superreads, region, **kwargs)
        )
    return all_descriptions


def assemble_io(
        input_superreads, input_resolution, output_assembly, minimum_weight=5,
        max_qs=2
):
    "IO function for performing assembly."
    superreads = read_json(input_superreads)
    resolution = read_json(input_resolution)
    assembly = assemble(
        superreads, resolution,
        minimum_weight=minimum_weight, max_qs=max_qs
    )
    write_json(output_assembly, assembly)


def local_reconstruction(
        superreads, assembly, consensus, covarying_sites
):
    "Reconstruct an individual region."
    number_of_quasispecies = np.sum(assembly['full_coverage'])
    number_of_covarying_sites = assembly['number_of_covarying_sites']
    covariation = [
        number_of_covarying_sites * ['N']
        for _ in range(number_of_quasispecies)
    ]
    counts = np.zeros((number_of_quasispecies, number_of_covarying_sites, 5))
    for qs_index, assembled in enumerate(assembly['original_indices']):
        for superread_index in assembled:
            superread = superreads[superread_index]
            weight = superread['weight']
            for vac_index, vac in enumerate(superread['vacs']):
                character_index = character_to_index_map[vac]
                site_index = superread['cv_start'] + vac_index
                counts[qs_index, site_index, character_index] += weight
        for site_index in np.arange(number_of_covarying_sites):
            character_index = np.argmax(counts[qs_index, site_index])
            character = index_to_character_map[character_index]
            covariation[qs_index][site_index] = character
    compressed_counts = np.sum(counts, axis=2)
    total_counts = np.sum(compressed_counts, axis=0)
    frequencies = pd.DataFrame((compressed_counts / total_counts).T)

    frequencies = np.sum(frequencies, axis=0) / frequencies.shape[0]
    frequencies = frequencies / np.sum(frequencies)
    quasispecies = []
    for i, qs_covariation in enumerate(covariation):
        label_content = 'quasispecies-%d_frequency-%.2f'
        label_data = (i, frequencies[i])
        label = label_content % label_data
        vacs = np.array(qs_covariation, dtype='<U1')
        seq = np.copy(consensus)
        seq[covarying_sites] = vacs
        sequence = Seq(''.join(seq))
        record = SeqRecord(
            seq=sequence, id=label, description=''
        )
        quasispecies.append(record)
    return quasispecies


def local_reconstruction_io(
        input_superreads, input_assembly, input_consensus,
        input_covarying_sites, region, output_quasispecies
    ):
    "IO function for local reconstruction."
    superreads = read_json(input_superreads)
    assembly = read_json(input_assembly)
    consensus = SeqIO.read(input_consensus, 'fasta')
    with open(input_covarying_sites) as json_file:
        covarying_sites = np.array(json.load(json_file), dtype=np.int)

    quasispecies = local_reconstruction(
        superreads, assembly[region], consensus, covarying_sites
    )
    SeqIO.write(quasispecies, output_quasispecies, 'fasta')


def ar_rate_estimation(superreads, description):
    "Estimate AR rate from leftover superreads."
    total_weight = sum([superread['weight'] for superread in superreads])
    ar_weight = 0
    ar_reads = []
    all_consensus = description['consensus']
    for superread in superreads:
        superread_length = superread['cv_end'] - superread['cv_start']
        if superread_length != len(superread['vacs']):
            print('WARNING! Discordance in superread data for:',
                  superread['index'])
            continue
        matches = [0 for _ in range(len(all_consensus))]
        for qs_index, qs_consensus in enumerate(all_consensus):
            cvs = range(superread['cv_start'], superread['cv_end'])
            for vacs_index, cvs_index in enumerate(cvs):
                qs_character = qs_consensus[cvs_index]
                superread_character = superread['vacs'][vacs_index]
                characters_agree = qs_character == superread_character
                if characters_agree:
                    matches[qs_index] += 1
        top_matches = sorted(matches, reverse=True)[:2]
        tag_as_ar = top_matches[0] > 0 and top_matches[1] > 0
        if tag_as_ar:
            ar_weight += superread['weight']
            ar_reads.append(superread['index'])
    ar_rate_estimate = ar_weight / total_weight
    return {
        'ar_rate_estimate': ar_rate_estimate,
        'ar_reads': ar_reads
    }


def ar_rate_estimation_io(input_superreads, input_description, output_json):
    "IO function for AR rate estimation."
    superreads = read_json(input_superreads)
    description = read_json(input_description)
    ar_info = ar_rate_estimation(superreads, description)
    write_json(output_json, ar_info)
