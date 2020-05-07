import json
import itertools as it

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy.linalg import lu_factor
from scipy.linalg import lu_solve


def simplex_projection(x):
    n = len(x)
    u = np.sort(x)[::-1]
    v = (np.cumsum(u)-1)/(np.arange(n)+1)
    K = np.argmax(v >= u)-1
    tau = v[K]
    return np.maximum(x-tau, np.zeros(n))


class QuadraticProx:
    def __init__(self, L, y, lambd):
        self.lambd = lambd
        LT = L.transpose()
        self.LTy = np.dot(LT, y)
        identity = np.eye(L.shape[1])
        A = identity + lambd*np.dot(LT, L)
        self.lu_factors = lu_factor(A)

    def apply(self, x):
        z = x + self.lambd*self.LTy
        return lu_solve(self.lu_factors, z)


def build_L_and_y(superreads, candidate_info):
    surviving_superread_indices = list(set(
        it.chain.from_iterable(candidate_info)
    ))
    surviving_superread_indices.sort()
    reindex_map = {ssi: i for i, ssi in enumerate(surviving_superread_indices)}
    surviving_superreads = [
        superreads[surviving_index]
        for surviving_index in surviving_superread_indices
    ]
    n = len(surviving_superread_indices)
    m = len(candidate_info)

    L = np.zeros((n, m))
    for i, row_info in enumerate(candidate_info):
        reindexed_row = [reindex_map[j] for j in row_info]
        L[reindexed_row, i] = 1

    y = np.array([
        row_info['frequency'] for row_info in surviving_superreads
    ])

    return L, y


def perform_regression(
        superread_info, candidate_info, maxit=100000, tolerance=1e-6, lambd=.1
        ):
    L, y = build_L_and_y(superread_info, candidate_info)
    quadratic_prox = QuadraticProx(L, y, lambd)

    y = np.ones(L.shape[1])/L.shape[1]
    x = simplex_projection(y)
    y = y + lambd*(quadratic_prox.apply(2*x - y) - x)
    for i in range(maxit):
        xold = np.array(x)
        x = simplex_projection(y)
        y = y + lambd*(quadratic_prox.apply(2*x - y) - x)
        relative_error = np.linalg.norm(x-xold)/np.linalg.norm(x)
        if i % 1000 == 0:
            params = (i, relative_error)
            print('Iteration %d, relative error %.2e...' % params)
        if relative_error < tolerance:
            params = (relative_error, i)
            print('Reached relative error of %.5e at iteration %d!' % params)
            break
    quasispecies_info = []
    for i, frequency in enumerate(x):
        if frequency > 0:
            print('Candidate %2d: frequency %.3f' % (i, frequency))
            quasispecies_info.append({
                'index': i,
                'frequency': frequency,
                'describing_superreads': candidate_info[i]
            })
    return quasispecies_info


def obtain_quasispecies(all_quasispecies_info, superreads, consensus, covarying_sites):
    all_quasispecies = []
    consensus_np = np.array(list(str(consensus.seq)), dtype='<U1')
    for i, quasispecies_info in enumerate(all_quasispecies_info):
        current_sequence = np.copy(consensus_np)
        for superread_index in quasispecies_info['describing_superreads']:
            cv_start = superreads[superread_index]['cv_start']
            vacs = superreads[superread_index]['vacs']
            n_sites = len(vacs)
            cv_end = cv_start + n_sites
            reference_indices = covarying_sites[cv_start: cv_end]
            current_sequence[reference_indices] = list(vacs)
        record_info = (i + 1, quasispecies_info['frequency'])
        record = SeqRecord(
            Seq(''.join(current_sequence)),
            id='quasispecies-%d_frequency-%.5f' % record_info,
            description=''
        )
        all_quasispecies.append(record)
    return all_quasispecies


def regression_io(
        input_superreads, input_describing_json,
        input_consensus, input_covarying_sites, output_fasta
        ):
    with open(input_superreads) as superread_file:
        superreads = json.load(superread_file)
    with open(input_describing_json) as candidates_file:
        describing = json.load(candidates_file)
    with open(input_covarying_sites) as json_file:
        covarying_sites = np.array(json.load(json_file), dtype=np.int)

    quasispecies_info = perform_regression(superreads, describing)
    consensus = SeqIO.read(input_consensus, 'fasta')
    all_quasispecies = obtain_quasispecies(
        quasispecies_info, superreads, consensus, covarying_sites
    )
    SeqIO.write(all_quasispecies, output_fasta, 'fasta')
