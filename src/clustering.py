import numpy as np
from sklearn.cluster import AgglomerativeClustering as AC
import matplotlib.pyplot as plt
import pickle

from skbio import DistanceMatrix
from skbio.tree import nj, TreeNode

from scipy.optimize import least_squares  # root

N = 132  # 270  # 66  # 294

with open('austr.pp.1.asd', 'rb') as f: # temp_asd.asd
    delta = pickle.load(f)

with open('temp_eig.data', 'rb') as f:
    lambdas, vecs = pickle.load(f)

with open('../test/haplotype_labels.txt') as f:
    lines = f.readlines()

labels = [l.split()[0] for l in lines]  # + ['WCD'] * 7

labels_d = [l.split()[1][:3] for l in lines]

for (i, l) in enumerate(labels.copy()):
    if l == 'ABO':
        labels[i] = labels_d[i]

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
        }

######
labels0 = np.array([l.split()[0] for l in lines])

where = np.where(np.logical_or(labels0 == 'ABO', labels0 == 'WCD'))[0]

labels = np.array(labels)[where]

######


colors = [colormap[l].lower() for l in labels]

arr = np.hstack((lambdas.reshape((N, 1)), vecs.T)).copy()
arr = sorted(arr, key=lambda x: x[0], reverse=True)
for i, v in enumerate(arr.copy()):
    arr[i] = np.sqrt(v[0])*v[1:]


npop = 5
offset = 2

print(len(arr))
arr = arr[:npop+offset]
arr = np.array(arr)
print(arr.shape)
arr = arr.T

clusterer = AC(n_clusters=npop, compute_full_tree=True)
labs = clusterer.fit_predict(arr)
labels = [hex(l)[-1].upper() for l in labs]

for p, q in [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]:
    fig, ax = plt.subplots()
    ax.scatter(arr.T[p], arr.T[q], c=colors, s=100)
    for i, txt in enumerate(labels):
        ax.annotate(txt, (arr.T[p, i], arr.T[q, i]))
    fig.savefig(f'whatever{p}{q}.svg')

plt.show()

ds = np.zeros((npop, npop))
coords = np.zeros((npop, npop+offset))
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
plt.pcolor(ds)
plt.show()

ids = list(map(str, range(npop)))
dm = DistanceMatrix(ds, ids)

tree = nj(dm)

outgroup = '2'

rt = TreeNode(name='rt')
temp = tree.find(outgroup).parent
temp.append(rt)
rt.append(tree.find(outgroup))
new_tree = tree.root_at('rt')

print(new_tree.ascii_art())

d_ind = np.zeros((npop, npop), dtype='int16')
constraints = []

def add_indices(tr: TreeNode, current=0):
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

def make_D(Dv):
    D = np.zeros((npop, npop))
    for i in range(npop):
        for j in range(npop):
            if i != j:
                D[i, j] = Dv[d_ind[i, j]]
    return D

print('------------')
print(ns)
Dv = [4, 3, 2, 1]
D = make_D(Dv)
print(D)

def make_b(D):
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
 
b = make_b(D)
print('------------')
print(b)
print(np.linalg.det(b))

T = 6.07**2/2

inits = np.zeros((npop-1,))
maxs = np.zeros((npop-1,))
mins = np.zeros((npop-1,))

for k in range(npop - 1):
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

print(inits)
print(mins)
print(maxs)

ls = sorted(lambdas, reverse=True)
ls = np.array(ls[:npop-1])

def dev(dv):
    D = make_D(dv)
    b = make_b(D)
    vals, vecs = np.linalg.eig(b)
    real_vals = 2*(1 - vals)/T**2
    real_vals = np.array(sorted(real_vals, reverse=True))[:npop-1]
    return real_vals - ls

res = least_squares(dev, inits, bounds=(mins, maxs), gtol=1e-10)
print(res)
