import json
import itertools as it

import numpy as np
import networkx as nx
import pandas as pd


def check_compatability(superread_i, superread_j, minimum_overlap=2):
    if not 'index' in superread_i or not 'index' in superread_j:
        return (False, 0)
    if superread_i['index'] == superread_j['index']:
        return (False, 0)
    i_cv_start = superread_i['cv_start']
    i_cv_end = superread_i['cv_end']
    j_cv_start = superread_j['cv_start']
    j_cv_end = superread_j['cv_end']
    start_before_start = i_cv_start <= j_cv_start
    start_before_end = j_cv_start < i_cv_end
    end_before_end = i_cv_end <= j_cv_end
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


def get_full_edge_list(G):
    edge_list = []
    for superread_index_i in G.nodes:
        superread_i = G.nodes[superread_index_i]
        for superread_index_j in G.nodes:
            superread_j = G.nodes[superread_index_j]
            should_include_edge, edge_data = check_compatability(
                superread_i, superread_j
            )
            if should_include_edge:
                edge_list.append({
                    'i': superread_i['index'],
                    'j': superread_j['index'],
                    **edge_data
                })
    return edge_list


def initialize_superread_graph(superreads, weight_percentile_cutoff=.1,
        minimum_weight=3, minimum_cvs=3
        ):
    G = nx.DiGraph()
    n_cv = max([sr['cv_end'] for sr in superreads])
    G.add_node('source', **{'cv_start': 0, 'cv_end': 0})
    G.add_node('target', **{'cv_start': n_cv, 'cv_end': n_cv})
    superread_weights = pd.Series([
        sr['weight']
        for sr in superreads
        if sr['weight'] >= minimum_weight
    ])
    weight_cutoff = superread_weights.quantile(weight_percentile_cutoff)
    message = 'Building superread graph with %d out of %d superreads...'
    admissible_count = (superread_weights >= weight_cutoff).sum()
    print(message % (admissible_count, len(superreads)))
    for superread in superreads:
        meets_cutoff = superread['weight'] >= weight_cutoff
        above_minimum = superread['weight'] >= minimum_weight
        enough_covariation = len(superread['vacs']) >= minimum_cvs
        admissible = meets_cutoff and above_minimum and enough_covariation
        if admissible:
            G.add_node(superread['index'], **superread)
            if superread['cv_start'] == 0:
                G.add_edge('source', superread['index'])
            if superread['cv_end'] == n_cv:
                G.add_edge(superread['index'], 'target')
    return G


def transitive_reduction(G):
    Gtr = nx.algorithms.dag.transitive_reduction(G)
    for node in Gtr.nodes:
        Gtr.nodes[node].update(G.nodes[node])
    for edge in Gtr.edges:
        source, target = edge
        Gtr.edges[source, target].update(G.edges[source, target])
    return Gtr


def create_full(superreads, weight_percentile_cutoff=.1,
        minimum_weight=3
    ):
    G = initialize_superread_graph(
        superreads, weight_percentile_cutoff, minimum_weight
    )
    all_edges = get_full_edge_list(G)
    for edge in all_edges:
        G.add_edge(edge['i'], edge['j'], weight=edge['overlap'])
    return transitive_reduction(G)


def dynamic_programming_path_count(G, source='source', target='target'):
  G.nodes[source]['npath'] = 1
  for node in nx.dfs_postorder_nodes(G, source):
    if node == 'target':
      G.nodes[node]['npath'] = 1
    else:
      number_of_current_paths = sum([
        G.nodes[successor]['npath']
        for successor in G.successors(node)
      ])
      G.nodes[node]['npath'] = number_of_current_paths
  return G.nodes[source]['npath']


def annotate_and_serialize_graph(G):
    superread_json = nx.node_link_data(G)
    superread_json['number_of_paths'] = dynamic_programming_path_count(G)
    superread_json['number_of_nodes'] = len(G)
    superread_json['number_of_edges'] = G.number_of_edges()
    return superread_json


def graph_io(input_srdata, output_json, graph_type='full',
        edges_per_node=3, weight_percentile_cutoff=.5,
        minimum_weight=5, edge_key='overlap'):
    with open(input_srdata) as json_file:
        superreads = json.load(json_file)
    if graph_type == 'full':
        G = create_full(superreads)
    elif graph_type == 'reduced':
        G = create_simple_reduced(
            superreads, edges_per_node, weight_percentile_cutoff,
            minimum_weight, edge_key
        )
    elif graph_type == 'incremental':
        G = create_incremental_reduced(
            superreads, weight_percentile_cutoff,
            minimum_weight, edge_key
        )
    else:
        raise ValueError('Encountered unknown graph type:i %s' % graph_type)
    superread_json = annotate_and_serialize_graph(G)
    with open(output_json, 'w') as json_file:
        json.dump(superread_json, json_file, indent=2)


def create_simple_reduced(superreads, edges_per_node=2,
        weight_percentile_cutoff=.1, minimum_weight=5, edge_key='overlap'
    ):
    G = initialize_superread_graph(
        superreads, weight_percentile_cutoff, minimum_weight
    )
    all_edges = get_full_edge_list(G)
    grouped_edges = it.groupby(all_edges, lambda edge: edge['i'])
    for i, edges_i in grouped_edges:
        sorted_edges = sorted(
            edges_i, reverse=True, key=lambda edge: edge[edge_key]
        )
        for edge in it.islice(sorted_edges, 0, edges_per_node):
            G.add_edge(edge['i'], edge['j'], weight=edge[edge_key])
    G_tr = transitive_reduction(G)
    number_of_paths = dynamic_programming_path_count(G_tr)
    message = 'Built superread graph with %d nodes, %d edges, and %d paths.'
    print(message % (len(G_tr.nodes), len(G_tr.edges), number_of_paths))
    return G_tr


def create_incremental_reduced(superreads, weight_percentile_cutoff=.1,
        minimum_weight=3, edge_key='overlap', edges_per_iteration=100, max_paths=100):
    G = initialize_superread_graph(
        superreads, weight_percentile_cutoff, minimum_weight
    )
    all_edges = sorted(
        get_full_edge_list(G),
        key=lambda edge: edge[edge_key],
        reverse=True
    )
    print('Building read graph incrementally...')
    for i, edge in enumerate(all_edges):
        G.add_edge(edge['i'], edge['j'], weight=edge[edge_key])
        if i % edges_per_iteration == 0:
            G = transitive_reduction(G)
            number_of_paths = dynamic_programming_path_count(G)
            print('Iteration:', i, ' -- Edges: ', len(G.edges), ' -- Paths: ', number_of_paths)
            if number_of_paths > max_paths:
                break
    number_of_paths = dynamic_programming_path_count(G)
    message = 'Built superread graph with %d nodes, %d edges, and %d paths.'
    print(message % (len(G.nodes), len(G.edges), number_of_paths))
    return G


def get_candidates(G):
    number_of_covarying_sites = G.nodes['target']['cv_end']
    G.nodes['source']['candidate_quasispecies'] = np.array(
        [[]], dtype='<U1'
    )
    G.nodes['source']['describing_superreads'] = [[]]
    G.nodes['target']['vacs'] = ''
    reverse_post_order = list(
        nx.dfs_postorder_nodes(G, 'source')
        )[::-1][1:]

    G.nodes['source']['npath'] = 1
    for node in reverse_post_order:
        for predecessor in G.predecessors(node):
            if not 'npath' in G.nodes[predecessor]:
                G.nodes[predecessor]['npath'] = 0
        number_of_current_paths = sum([
            G.nodes[predecessor]['npath']
            for predecessor in G.predecessors(node)
        ])
        G.nodes[node]['npath'] = number_of_current_paths
    number_of_candidates = G.nodes['target']['npath']
    if number_of_candidates > 10000:
        error_message = '%d candidates... refusing to proceed!'
        raise ValueError(error_message % number_of_candidates)
    print(
        'Obtaining',
        number_of_candidates,
        'candidate quasispecies...'
        )

    for descendant in reverse_post_order:
        descendant_start = G.nodes[descendant]['cv_start']
        vacs = G.nodes[descendant]['vacs']
        extended_candidates = []
        extended_describing_superreads = []
        for ancestor in G.pred[descendant].keys():
            if not 'candidate_quasispecies' in G.nodes[ancestor]:
                continue
            ancestor_end = G.nodes[ancestor]['cv_end']
            vacs_start_index = ancestor_end - descendant_start
            vacs_np = np.array(
                [list(vacs[vacs_start_index:])], dtype='<U1'
                )
            ancestor_candidates = \
                G.nodes[ancestor]['candidate_quasispecies']
            n_ancestor_candidates = ancestor_candidates.shape[0]
            current_extended_candidates = np.hstack([
                ancestor_candidates,
                np.repeat(vacs_np, n_ancestor_candidates, axis=0)
            ])
            extension = [descendant] if descendant != 'target' else []
            current_extended_describing_superreads = [
                describing_list + extension
                for describing_list
                in G.nodes[ancestor]['describing_superreads']
            ]
            extended_candidates.append(current_extended_candidates)
            extended_describing_superreads.append(
                current_extended_describing_superreads
            )
        candidate_quasispecies = np.vstack(extended_candidates)
        G.nodes[descendant]['candidate_quasispecies'] = \
            candidate_quasispecies
        G.nodes[descendant]['describing_superreads'] = list(
            it.chain.from_iterable(extended_describing_superreads)
        )
    candidate_quasispecies = G.nodes['target']['candidate_quasispecies']
    describing_superreads = G.nodes['target']['describing_superreads']
    return candidate_quasispecies, describing_superreads


def deserialize_graph(input_graph_json):
    with open(input_graph_json) as json_file:
        graph_json = json.load(json_file)
    G = nx.node_link_graph(graph_json)
    return G


def candidates_io(input_graph_json, output_describing_json):
    G = deserialize_graph(input_graph_json)
    candidate_vacs, describing_superreads = get_candidates(G)
    with open(output_describing_json, 'w') as json_file:
        json.dump(describing_superreads, json_file, indent=2)
