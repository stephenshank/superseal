import json
import itertools as it

import networkx as nx


def check_compatability(superread_i, superread_j, minimum_overlap=2):
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
        long_enough = overlap >= minimum_overlap
        compatible = agree_on_overlap and long_enough
        return (compatible, overlap)
    return (False, 0)


def create_full(superreads, covarying_sites):
    G = nx.DiGraph()
    G.add_node('source')
    G.add_node('target')
    n_cv = len(covarying_sites)
    for superread in superreads:
        G.add_node(superread['index'], **superread)
        if superread['cv_start'] == 0:
            G.add_edge('source', superread['index'])
        if superread['cv_end'] == n_cv:
            G.add_edge(superread['index'], 'target')
    for superread_i in superreads:
        for superread_j in superreads:
            should_include_edge, overlap = check_compatability(
                superread_i, superread_j
            )
            if should_include_edge:
                G.add_edge(
                    superread_i['index'], superread_j['index'],
                    overlap=overlap
                )
    return G


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


def full_graph_io(input_srdata, input_cvs, output_json):
    with open(input_cvs) as json_file:
        covarying_sites = json.load(json_file)
    with open(input_srdata) as json_file:
        superreads = json.load(json_file)
    G = create_full(superreads, covarying_sites)
    G = nx.algorithms.dag.transitive_reduction(G)
    superread_json = nx.node_link_data(G)
    superread_json['number_of_paths'] = dynamic_programming_path_count(G)
    with open(output_json, 'w') as json_file:
        json.dump(superread_json, json_file, indent=2)
