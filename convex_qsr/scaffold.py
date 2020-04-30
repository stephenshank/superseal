import json

import numpy as np
import pandas as pd

from .graph import check_compatability


class Scaffold:
    def __init__(self, superread_i, superread_j, number_of_covarying_sites):
        self.members = {}
        self.coverage = np.zeros(number_of_covarying_sites, dtype=np.int)
        self.merge_node(superread_i)
        self.merge_node(superread_j)

    def merge_node(self, superread):
        indices = np.arange(superread['cv_start'], superread['cv_end'])
        self.coverage[indices] += superread['weight']
        self.members[superread['filtered_index']] = True

    def handle_coverage(self, superread_i, superread_j):
        self.handle_single_coverage(superread_i)
        self.handle_single_coverage(superread_j)

    def merge_edge(self, superread_i, superread_j):
        self.merge_node(superread_i)
        self.merge_node(superread_j)

    def merge_scaffold(self, scaffold):
        for i in scaffold.members.keys():
            self.members[i] = True
        self.coverage += scaffold.coverage

    def check_membership(self, superread):
        return superread['filtered_index'] in self.members

    def check_membership_by_index(self, superread_index):
        return superread_index in self.members


def filter_superreads(superreads, minimum_weight):
    filtered_superreads = list(filter(
        lambda superread: superread['weight'] >= minimum_weight,
        superreads
    ))
    number_of_covarying_sites = 0
    for i, superread in enumerate(filtered_superreads):
        superread['filtered_index'] = i
        number_of_covarying_sites = max(
            number_of_covarying_sites, superread['cv_end']
        )
    return filtered_superreads, number_of_covarying_sites


def get_edge_list(superreads):
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
                    **edge_data
                })
    return pd.DataFrame(edge_list).sort_values(by='support', ascending=False)


def check_integrity(scaffolds):
    for i, scaffold_i in enumerate(scaffolds):
        for member in scaffold_i.members.keys():
            for scaffold_j in scaffolds[i+1:]:
                assert not scaffold_j.check_membership_by_index(member)


def scaffold_qsr(all_superreads, minimum_weight=5, max_qs=2, verbose=True):
    superreads, number_of_covarying_sites = filter_superreads(
        all_superreads, minimum_weight
    )
    edges = get_edge_list(superreads)
    scaffolds = []
    current_edge = 0
    for _, edge in edges.iterrows():
        if verbose:
            print('EDGE %d OF %d...' % (current_edge, len(edges)))
            current_edge += 1
        superread_i = superreads[edge['i']]
        superread_j = superreads[edge['j']]
        i_in = None
        j_in = None
        for scaffold_index, scaffold in enumerate(scaffolds):
            if scaffold.check_membership(superread_i):
                if not i_in is None:
                    Exception('Superread found in distinct scaffolds!')
                i_in = scaffold_index
            if scaffold.check_membership(superread_j):
                if not j_in is None:
                    Exception('Superread found in distinct scaffolds!')
                j_in = scaffold_index
        completely_new = i_in is None and j_in is None
        only_i_connects = not i_in is None and j_in is None
        only_j_connects = i_in is None and not j_in is None
        i_and_j_connect = not i_in is None and not j_in is None
        i_and_j_in_same = i_and_j_connect and i_in == j_in
        i_and_j_in_distinct = i_and_j_connect and i_in != j_in
        if completely_new:
            message = 'Creating new scaffold from superreads %d and %d...'
            message_data = (edge['i'], edge['j'])
            new_scaffold = Scaffold(
                superread_i, superread_j, number_of_covarying_sites
            )
            scaffolds.append(new_scaffold)
        elif only_i_connects or i_and_j_in_same:
            message = 'Merging superreads %d and %d to scaffold %d...'
            message_data = (edge['i'], edge['j'], i_in)
            scaffolds[i_in].merge_edge(superread_i, superread_j)
        elif only_j_connects:
            message = 'Merging superreads %d and %d to scaffold %d...'
            message_data = (edge['i'], edge['j'], j_in)
            scaffolds[j_in].merge_edge(superread_i, superread_j)
        elif i_and_j_in_distinct:
            message = 'Merging scaffolds %d and %d, and edge %d and %d...'
            smaller_index = min(i_in, j_in)
            larger_index = max(i_in, j_in)
            message_data = (smaller_index, larger_index, edge['i'], edge['j'])
            scaffolds[smaller_index].merge_edge(superread_i, superread_j)
            scaffolds[smaller_index].merge_scaffold(scaffolds[larger_index])
            del scaffolds[larger_index]
        if verbose:
            print(message % message_data)
        fully_covered = 0
        for scaffold in scaffolds:
            if np.all(scaffold.coverage > 0):
                fully_covered += 1
        if verbose:
            message = '%d total scaffolds, %d which fully cover...'
            message_data = (len(scaffolds), fully_covered)
            print(message % message_data)
        if fully_covered >= max_qs:
            print(fully_covered, '>=', max_qs, 'quasispecies... terminating!')
            break
        check_integrity(scaffolds)
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
        ]
    }


def scaffold_qsr_io(input_superreads, output_describing):
    with open(input_superreads) as json_file:
        superreads = json.load(json_file)
    describing_superreads = scaffold_qsr(superreads)
    with open(output_describing, 'w') as json_file:
        json.dump(describing_superreads, json_file, indent=2)
