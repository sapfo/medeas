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


def get_coordinate(simulation, p):
    with open(simulation.vec_pattern.format(p), 'rb') as f:
        lambdas, vecs = pickle.load(f)
    N = len(lambdas)

    coordinates = np.hstack((lambdas.reshape((N, 1)), vecs.T)).copy()
    coordinates = sorted(coordinates, key=lambda x: x[0], reverse=True)

    for i, v in enumerate(coordinates.copy()):
        coordinates[i] = np.sqrt(v[0])*v[1:]

    npop = simulation.K
    coordinates = coordinates[:npop + OFFSET]
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

    clusterer = AC(n_clusters=npop, compute_full_tree=True)
    lab_infered = clusterer.fit_predict(coordinates)

    return lab_infered


def find_tree(npop: int, asd_file: str,
              inferred_labels: 'np.ndarray[int]',
              arr: 'np.ndarray[float]',
              simulation,
              ) -> Tuple[TreeNode, 'np.ndarray[int]', 'np.ndarray[float]']:
    """Find tree topology using the centers of mass of clusters.
    'inferred_labels' contains assigned labels. Return the neighbor join tree, population sizes,
    and the bloks of original distance matrix that correspond to given
    population pairs (for further determination of fitting window).
    """

    with open(asd_file, 'rb') as f:
        delta = pickle.load(f)

    ds = np.zeros((npop, npop))
    coords = np.zeros((npop, npop+OFFSET))
    ns = np.zeros((npop,))

    for i in set(inferred_labels):
        coords[i, :] = np.mean(arr[np.where(inferred_labels == i)[0], :], axis=0)
        ns[i] = len(np.where(inferred_labels == i)[0])

    blocks = np.zeros((npop, npop), dtype='object')

    for i in range(npop):
        for j in range(npop):
            blocks[i, j] = delta[np.where(inferred_labels == i)[0]].T[np.where(inferred_labels == j)[0]]
    if simulation.output_level >= 1:
        print(coords)
        print(coords.shape)

    for i in range(npop):
        for j in range(npop):
            ds[i, j] = np.sqrt(np.sum((coords[i] - coords[j])**2))
    if simulation.output_level >= 1:
        print(ds)
    if simulation.output_level:
        plt.pcolor(ds)
        plt.show()

    ids = list(map(str, range(npop)))
    dm = DistanceMatrix(ds, ids)
    if npop == 2:
        tree = read(StringIO('(0:0.1, 1:0.1);'), format='newick', into=TreeNode)
        return tree, ns, blocks
    tree = nj(dm)
    new_tree = tree.root_at_midpoint()
    print(new_tree.ascii_art())
    return new_tree, ns, blocks


def find_distances(npop: int, T: float,
                   new_tree: TreeNode, ns: 'np.ndarray[int]',
                   lambdas: 'np.ndarray[float]',
                   blocks: 'np.ndarray[np.ndarray[float]]',
                   simulation
                   ) -> Tuple[OptimizeResult, List[Tuple[int, int]]]:
    """Find split times from the tree topology."""
    d_ind = -np.ones((npop, npop), dtype='int16')
    constraints = []

    def add_indices(tr: TreeNode, current: int = 0) -> None:
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
                d_ind[i, j] = d_ind[j, i] = current
        previous = current
        current += 1
        if not left.is_tip():
            constraints.append((current, previous))
            add_indices(left, current)
            current += 1 + sum([1 for _ in left.non_tips()])
        if not right.is_tip():
            constraints.append((current, previous))
            add_indices(right, current)

    add_indices(new_tree)

    if simulation.output_level == 2:
        print(d_ind)
        print(constraints)

    def make_D(Dv: List[int]) -> 'np.ndarray[int]':
        """Make distance (split time) matrix from the given vector of
        split times 'Dv'.
        """
        D = np.zeros((npop, npop))
        for i in range(npop):
            for j in range(npop):
                if i != j:
                    D[i, j] = Dv[d_ind[i, j]]
        return D

    if simulation.output_level == 2:
        print('------------')
        print(ns)
        Dv = range(npop, 0, -1)
        D = make_D(Dv)
        print(D)

    def make_b(D: 'np.ndarray[int]') -> 'np.ndarray[flost]':
        """Make the B-matrix (see paper for details)."""
        n = np.sum(ns)
        b = np.zeros((npop, npop))
        delta = np.zeros((npop,))
        for i in range(npop):
            delta[i] = (sum(ns[k]*(D[i, k] + 1)**2 for k in range(npop) if k != i) + ns[i] - 1)/n
        delta_0 = sum(ns[i]*delta[i] for i in range(npop))/n
        for i in range(npop):
            for j in range(npop):
                b[i, j] = ns[i] * ((D[i, j] + 1)**2 - delta[i] - delta[j] + delta_0)
        return b

    if simulation.output_level == 2:
        b = make_b(D)
        print('------------')
        print(b)
        print(np.linalg.det(b))

    T = T**2/2

    inits = np.zeros((npop-1,))
    maxs = np.zeros((npop-1,))
    mins = np.zeros((npop-1,))

    for k in range(npop - 1):
        # From every block of distances find boundaries for fitting,
        # i.e. smallest and largest value.
        ks = np.where(d_ind == k)
        subblocks = blocks[ks]
        means = []
        s_mins = []
        s_maxs = []
        for sb in subblocks:
            means.append(np.mean(sb))
            s_mins.append(sb.min())
            s_maxs.append(sb.max())
        inits[k] = np.mean(means)
        mins[k] = min(s_mins)
        maxs[k] = max(s_maxs)

    inits = T*inits/2 - 1
    print(inits)
    mins = T*mins/2 - 1
    maxs = T*maxs/2 - 1
    for i, (mn, mx) in enumerate(zip(mins, maxs)):
        inits[i] = uniform(mn, mx)

    if simulation.output_level == 2:
       print('Inits:', inits)
       print('Mins:', mins)
       print('Maxs:', maxs)

    ls = sorted(lambdas, reverse=True)
    ls = np.array(ls[:npop-1])

    def dev(dv: 'np.ndarray[float]') -> 'np.ndarray[float]':
        """Find deviation (residuals) for given split time vector 'dv'."""
        D = make_D(dv)
        b = make_b(D)
        vals, vecs = np.linalg.eigh(b)
        real_vals = 2*(1 - vals)/T**2
        real_vals = np.array(sorted(real_vals, reverse=True))[:npop-1]
        #Even if we expect the eigenvalue to be real, the "numerical" tiny little complex part cause problem for some simulations
        return np.real(real_vals - ls)

    res = least_squares(dev, inits, bounds=(mins, maxs), gtol=1e-15)
    return res, constraints


def validate_dists(dists: 'np.ndarray[float]',
                   constraints: List[Tuple[int, int]]) -> bool:
    """Verify if the given split times satisfy the 'constrains' obtained from
    the tree topology.
    """
    for c in constraints:
        if not smaller(dists[c[0]], dists[c[1]]):
            return False
    for dist in dists:
        if dist < -0.01:
            return False
    return True


def smaller(x: float, y: float) -> bool:
    """Approximate 'less than equal'."""
    if x - y < 0.01:
        return True
    if x/y < 1.01 and np.sign(y) > np.sign(x):
        return True
    return False