import numpy as np
from sklearn.cluster import AgglomerativeClustering as AC
import matplotlib.pyplot as plt
import pickle

from skbio import DistanceMatrix
from skbio.tree import nj, TreeNode
from io import StringIO
from skbio import read

from scipy.optimize import OptimizeResult, least_squares  # or root?

from typing import Tuple, List
from collections import Counter
from random import uniform

OFFSET = 2
from options import TESTING


def plot_mds(npop: int,
                       vectors_file: str, labels_file: str
                       ) -> Tuple['np.ndarray[int]', 'np.ndarray[float]',
                                  'np.ndarray[float]', List[str]]:
    """Perform agglomerative clustering for 'npop' clusters.

    'vectors_file' and 'labels_file' are names of files to read data.
    Return found labels, distaces from large eigenvalues,
    eigenvalues read from file, and labels read from file.
    """

    with open(labels_file) as f:
        lines = f.readlines()

    labels = [l.split()[0] for l in lines]  # + ['WCD'] * 7

    labels_d = [l.split()[1][:3] for l in lines]

    for (i, l) in enumerate(labels.copy()):
        if l == 'ABO':
            labels[i] = labels_d[i]  # TODO: maybe move this logic to main.py
                                     # or a separate function?

    colormap = {
    'CAI': '#2679B2',
    'WPA': '#CAB3D5',
    'BDV': '#A7CEE2',
    'RIV': '#E01F27',
    'WCD': '#FCBE75',
    'WON': '#FD7F23',
    'ENY': '#B3DE8E',
    'NGA': '#399F34',
    'PIL': '#F99B9B',
    'CHI': 'black',
    'BRI': 'blue',
    'PAP': 'red'
    }  # TODO: Autogenerate colormap

    colors = [colormap[l].lower() for l in labels]
    res_labels = labels.copy()

    arr = np.hstack((lambdas.reshape((N, 1)), vecs.T)).copy()
    arr = sorted(arr, key=lambda x: x[0], reverse=True)
    for i, v in enumerate(arr.copy()):
        arr[i] = np.sqrt(v[0])*v[1:]

    print(len(arr))
    arr = arr[:npop+OFFSET]
    arr = np.array(arr)
    arr = arr.T
    print('clustering will be performed on a ' + str(arr.shape) + ' matrix')

    clusterer = AC(n_clusters=npop, compute_full_tree=True)
    labs = clusterer.fit_predict(arr)
    labels = [hex(l)[-1].upper() for l in labs]

    if TESTING:
        # TODO: autogenerate nice summary plot depending on 'npop'
        for p, q in [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]:
            fig, ax = plt.subplots()
            ax.scatter(arr.T[p], arr.T[q], c=colors, s=100)
            for i, txt in enumerate(labels):
                ax.annotate(txt, (arr.T[p, i], arr.T[q, i]))
            fig.savefig(f'whatever{p}{q}.svg')

        plt.show()
    return labs, arr, lambdas, res_labels


def find_tree(npop: int, asd_file: str,
              labs: 'np.ndarray[int]',
              arr: 'np.ndarray[float]',
              outgroups: List[str], res_labels: List[str]
              ) -> Tuple[TreeNode, 'np.ndarray[int]', 'np.ndarray[float]']:
    """Find tree topology using the centers of mass of clusters.

    'labs' contains assigned labels. 'asd_file' is the name of the file to read
    original distance matrix. Return the neighbor join tree, population sizes,
    and the bloks of original distance matrix that correspond to given
    population pairs (for further determination of fitting window).
    """
    # TODO: Refactor functions in this file to be more logical instead of
    # carrying over unrelated info.

    res_labels = np.array(res_labels)
    cond_lab = np.zeros(labs.shape)
    for outg in outgroups:
        cond_lab = np.logical_or(cond_lab, res_labels == outg)
    outg_labs = labs[np.where(cond_lab)]
    count = Counter(outg_labs)
    if len(outg_labs):
        outgroup = count.most_common()[0][0]
        outgroup = hex(outgroup)[-1].upper()

    with open(asd_file, 'rb') as f:
        delta = pickle.load(f)

    ds = np.zeros((npop, npop))
    coords = np.zeros((npop, npop+OFFSET))
    ns = np.zeros((npop,))

    for i in set(labs):
        coords[i, :] = np.mean(arr[np.where(labs == i)[0], :], axis=0)
        ns[i] = len(np.where(labs == i)[0])

    blocks = np.zeros((npop, npop), dtype='object')

    for i in range(npop):
        for j in range(npop):
            blocks[i, j] = delta[np.where(labs == i)[0]].T[np.where(labs == j)[0]]

    print(coords)
    print(coords.shape)

    for i in range(npop):
        for j in range(npop):
            ds[i, j] = np.sqrt(np.sum((coords[i] - coords[j])**2))

    print(ds)
    if TESTING:
        plt.pcolor(ds)
        plt.show()

    ids = list(map(str, range(npop)))
    dm = DistanceMatrix(ds, ids)
    if npop == 2:
        tree = read(StringIO('(0:0.1, 1:0.1);'), format='newick', into=TreeNode)
        return tree, ns, blocks
    tree = nj(dm)

    rt = TreeNode(name='rt')
    temp = tree.find(outgroup).parent
    temp.append(rt)
    rt.append(tree.find(outgroup))
    new_tree = tree.root_at('rt')

    print(new_tree.ascii_art())
    return new_tree, ns, blocks


def find_distances(npop: int, T: float,
                   new_tree: TreeNode, ns: 'np.ndarray[int]',
                   lambdas: 'np.ndarray[float]',
                   blocks: 'np.ndarray[np.ndarray[float]]'
                   ) -> Tuple[OptimizeResult, List[Tuple[int, int]]]:
    """Find split times from the tree topology."""
    d_ind = np.zeros((npop, npop), dtype='int16')
    constraints = []

    def add_indices(tr: TreeNode, current: int = 0) -> None:
        """For every pair of populations substitute an index of coresponding split time.
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
                d_ind[i, j] = current
        previous = current
        current += 1
        if not left.is_tip():
            constraints.append((current, previous))
            add_indices(left, current)
            current += 1
        if not right.is_tip():
            constraints.append((current, previous))
            add_indices(right, current)

    add_indices(new_tree)
    d_ind = d_ind + d_ind.T
    for i in range(npop):
        d_ind[i, i] = -1
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

    if TESTING:
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

    if TESTING:
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
    mins = T*mins/2 - 1
    maxs = T*maxs/2 - 1
    for i, (mn, mx) in enumerate(zip(mins, maxs)):
        inits[i] = uniform(mn, mx)

    if TESTING:
        print('Inits:', inits)
        print('Mins:', mins)
        print('Maxs:', maxs)

    ls = sorted(lambdas, reverse=True)
    ls = np.array(ls[:npop-1])

    def dev(dv: 'np.ndarray[float]') -> 'np.ndarray[float]':
        """Find deviation (residuals) for given split time vector 'dv'."""
        D = make_D(dv)
        b = make_b(D)
        vals, vecs = np.linalg.eig(b)
        real_vals = 2*(1 - vals)/T**2
        real_vals = np.array(sorted(real_vals, reverse=True))[:npop-1]
        return real_vals - ls

    # TODO: find more solutions for D using different initial points (grid(-) or random(+))
    res = least_squares(dev, inits, bounds=(mins, maxs), gtol=1e-15)
    return res, constraints


# TODO: implement validations: positive values and order
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
    if y -  x < 0.001:
        return True
    if x/y < 1.001:
        return True
    return False