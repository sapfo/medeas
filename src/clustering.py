import numpy as np
import os
from sklearn.cluster import AgglomerativeClustering as AC
import matplotlib.pyplot as plt
import pickle

from skbio import DistanceMatrix
from skbio.tree import nj, TreeNode
from io import StringIO
from skbio import read

from scipy.optimize import OptimizeResult, least_squares  # or root?

from typing import Tuple, List
from random import uniform

OFFSET = 2


def get_mds_coordinate(simulation, p):
    """Compute the MDS coordinate from the eigenvalues and eigenvectore of the
    MDS matrix"""
    with open(simulation.vec_pattern.format(p), 'rb') as f:
        lambdas, vecs = pickle.load(f)
    N = len(lambdas)

    coordinates = np.hstack((lambdas.reshape((N, 1)), vecs.T)).copy()
    coordinates = sorted(coordinates, key=lambda x: x[0], reverse=True)

    for i, v in enumerate(coordinates.copy()):
        if v[0] > 0:
            coordinates[i] = np.sqrt(v[0])*v[1:]
        else:
            coordinates[i] = 0 * v[1:]


    coordinates = np.array(coordinates)
    coordinates = coordinates.T
    return(coordinates)

def perform_clustering(npop: int,
                       coordinates,
                       simulation
                       ) -> Tuple['np.ndarray[int]', 'np.ndarray[float]',
                                  'np.ndarray[float]']:
    """Perform agglomerative clustering for simulation
    Return found labels, distances from large eigenvalues,
    eigenvalues read from file, and labels read from file.
    """

    if simulation.output_level >= 1:
        print('clustering will be performed on a ' + str(coordinates.shape) + ' matrix')

    clusterer = AC(n_clusters=npop, compute_full_tree=True,linkage="ward")
    lab_infered = clusterer.fit_predict(coordinates)

    return lab_infered
def build_distance_subblock(npop, labels, delta):
    blocks = np.zeros((npop, npop), dtype='object')

    for i in range(npop):
        for j in range(npop):
            blocks[i, j] = delta[np.where(labels == i)[0]].T[np.where(labels == j)[0]]
    return(blocks)

def build_population_dimension(npop: int, numerical_labels):
    ns = np.zeros((npop,),dtype=int)
    for i in set(numerical_labels):
        ns[i] = len(np.where(numerical_labels == i)[0])
    return(ns)

def find_tree(npop: int,
              numerical_label: 'np.ndarray[int]',
              arr: 'np.ndarray[float]',
              ) -> TreeNode:
    """Find tree topology using the centers of mass of clusters.
    'inferred_labels' contains assigned labels. Return the neighbor join tree, population sizes,
    and the bloks of original distance matrix that correspond to given
    population pairs (for further determination of fitting window).
    """
    if npop == 2:
        tree = read(StringIO('(0:0.1, 1:0.1);'), format='newick', into=TreeNode)
        return tree

    arr = arr[:, :npop + OFFSET]
    ds = np.zeros((npop, npop))
    coords = np.zeros((npop, npop+OFFSET))
    for i in set(numerical_label):
        coords[i, :] = np.mean(arr[np.where(numerical_label == i)[0], :], axis=0)
    for i in range(npop):
        for j in range(npop):
            ds[i, j] = np.sqrt(np.sum((coords[i] - coords[j])**2))

    ids = list(map(str, range(npop)))
    dm = DistanceMatrix(ds, ids)
    tree = nj(dm)
    new_tree = tree.root_at_midpoint()
    print(new_tree.ascii_art())
    print(new_tree)
    return new_tree

def set_tree_from_input(asd_file, simulation) -> Tuple[TreeNode, 'np.ndarray[int]', 'np.ndarray[float]']:
    """Using the given tree topology, Return the neighbor join tree, population sizes,
    and the bloks of original distance matrix that correspond to given
    population pairs (for further determination of fitting window).
    """
    print(simulation.topology)
    tree = read(StringIO(simulation.topology),format='newick', into=TreeNode)
    print(tree.ascii_art())
    return tree

def build_split_index_matrix(tr: TreeNode, split_index_matrix: 'np.ndarray[int]', constraints: List, constraints_coal_time: List, current: int = 0) -> None:
    """For every pair of populations substitute an index of corresponding split time.
    Append distance constraints that correspond to given tree."""
    left = tr.children[0]
    right = tr.children[1]
    if left.is_tip():
        l_ind = [int(left.name)]
    else:
        l_ind = list(map(lambda x: int(x.name), left.tips()))
    if right.is_tip():
        r_ind = [int(right.name)]
    else:
        r_ind = list(map(lambda x: int(x.name), right.tips()))
    for i in l_ind:
        for j in r_ind:
            split_index_matrix[i, j] = split_index_matrix[j, i] = current
    previous = current
    current += 1
    if not left.is_tip():
        constraints.append((current, previous))
        build_split_index_matrix(left, split_index_matrix, constraints, constraints_coal_time, current)
        current += 1 + sum([1 for _ in left.non_tips()])
    else:
        constraints_coal_time.append((previous,left.name))
    if not right.is_tip():
        constraints.append((current, previous))
        build_split_index_matrix(right, split_index_matrix, constraints, constraints_coal_time, current)
    else:
        constraints_coal_time.append((previous, right.name))

def find_distances(npop: int, T: float, t_within: 'np.ndarray[float]',
                   tree: TreeNode, ns: 'np.ndarray[int]',
                   lambdas: 'np.ndarray[float]',
                   blocks: 'np.ndarray[np.ndarray[float]]',
                   output_level
                   ) -> Tuple[OptimizeResult, List[Tuple[int, int]]]:
    """Find coalescence times from the tree topology."""
    split_index_matrix = -np.ones((npop, npop), dtype='int16')
    constraints = []
    constraints_coal_time = []

    build_split_index_matrix(tree, split_index_matrix, constraints, constraints_coal_time)

    if output_level == 2:
        print(f"split_index_matrix = {split_index_matrix}")
        print(f"constraints_coal_time: {constraints_coal_time}")
    def make_tij(ts: List[int]) -> 'np.ndarray[float]':
        """Make distance (split time) matrix from the given vector of
        split times 'Dv'.
        """
        tij = np.ones((npop, npop))
        for i in range(npop):
            for j in range(npop):
                if i != j:
                    tij[i, j] = ts[split_index_matrix[i, j]]
                else:
                    tij[i, i] = t_within[i]
        return tij



    def make_b(ts: 'np.ndarray[int]') -> 'np.ndarray[float]':
        """Make the B-matrix (see paper for details)."""
        n = np.sum(ns)
        b = np.zeros((npop, npop))
        delta = np.zeros((npop,))
        for i in range(npop):
            delta[i] = (sum(ns[k]*ts[i, k]**2 for k in range(npop) if k != i) + ts[i,i]**2*(ns[i] - 1))/n
        delta_0 = sum(ns[i]*delta[i] for i in range(npop))/n
        for i in range(npop):
            for j in range(npop):
                if i != j:
                    b[i, j] = np.sqrt(ns[i]*ns[j]) * (ts[i, j]**2 - delta[i] - delta[j] + delta_0)
                else:
                    b[i, j] = (ns[i]-1) * ts[i, j] ** 2 - ns[i] * (2*delta[i] - delta_0)
        return b

    inits = np.zeros((npop-1,))
    mean = np.zeros((npop - 1,))
    maxs = np.zeros((npop-1,))
    mins = np.zeros((npop-1,))

    for k in range(npop - 1):
        # From every block of distances find boundaries for fitting,
        # i.e. smallest and largest value.
        ks = np.where(split_index_matrix == k)
        subblocks = blocks[ks]
        means = []
        s_mins = []
        s_maxs = []
        for sb in subblocks:
            means.append(np.mean(sb))
            s_mins.append(sb.min())
            s_maxs.append(sb.max())
        mean[k] = np.mean(means)
        mins[k] = min(s_mins)
        maxs[k] = max(s_maxs)

    mean = mean / 2
    mins = T*mins/2
    maxs = T*maxs/2
    for i, (minimum, maximum) in enumerate(zip(mins, maxs)):
        inits[i] = uniform(minimum, maximum)
    if output_level == 2:
        print(f"{100*'-'}")
        print(f"ns = {ns}")
        print('Inits:', inits)
        print('Mean:', mean)
        print('Mins:', mins)
        print('Maxs:', maxs)
        tij = make_tij(inits)
        print(f"{25*'-'} \n tij = \n {tij}")
        b = make_b(tij)
        print(f"{25*'-'} \n b = \n {b}")
        print(f"det (b) = {np.linalg.det(b)}")
        print(f"{100*'-'}")



    ls = sorted(lambdas, reverse=True)
    ls = np.array(ls[:npop-1])

    def dev(dv: 'np.ndarray[float]') -> 'np.ndarray[float]':
        """Find deviation (residuals) for given split time vector 'dv'."""
        D = make_tij(dv)
        b = make_b(D)
        vals, vecs = np.linalg.eigh(b)
        real_vals = -2*vals/T**2
        real_vals = np.array(sorted(real_vals, reverse=True))[:npop-1]
        return np.real(real_vals - ls)

    res = least_squares(dev, inits, bounds=(mins, maxs), gtol=1e-15)
    if output_level == 2:
        print('res:', res.x)
        tij = make_tij(res.x)
        print(f"{25*'-'} \n tij = \n {tij}")
        b = make_b(tij)
        print(f"{25*'-'} \n b = \n {b}")
        print(f"{100*'-'}")

    return res, constraints, constraints_coal_time


def validate_dists(dists: 'np.ndarray[float]',
                   ts: 'np.ndarray[float]',
                   constraints: List[Tuple[int, int]],
                   constraints_coal: List[Tuple[int, str]]
                   ) -> bool:
    """Verify if the given split times satisfy the 'constrains' obtained from
    the tree topology.
    """
    print(f"ts: {ts}")
    for c in constraints:
        if not smaller(dists[c[0]], dists[c[1]]):
            return False
    for c in constraints_coal:
        tts = ts[int(c[1])]
        d = dists[c[0]]
        if not smaller(tts, d):
        #if not smaller(dists[c[0]],ts[int(c[1])]):
            return False
    for dist in dists:
        if dist < -0.01:
            return False
    return True


def smaller(x: float, y: float) -> bool:
    """Approximate 'less than equal'."""
    if x - y < 0.000001:
        return True
    if x/y < 1.000001 and np.sign(y) > np.sign(x):
        return True
    return False