import numpy as np
from sklearn.cluster import AgglomerativeClustering as AC
from collections import defaultdict
import matplotlib.pyplot as plt
import pickle
from collections import Counter

N = 294

with open('temp_asd.asd', 'rb') as f:
    delta = pickle.load(f)

with open('temp_eig.data', 'rb') as f:
    lambdas, vecs = pickle.load(f)

arr = np.hstack((lambdas.reshape((N, 1)), vecs.T)).copy()
arr = sorted(arr, key=lambda x: x[0], reverse=True)
for i, v in enumerate(arr.copy()):
    arr[i] = np.sqrt(v[0])*v[1:]

print(len(arr))
arr = arr[:7]
arr = np.array(arr)
print(arr.shape)
arr = arr.T

npop = 8

clusterer = AC(n_clusters=npop)
labs = clusterer.fit_predict(arr)
labels = [hex(l)[-1].upper() for l in labs]

with open('../test/bla.tfam') as f:
    new_data = f.readlines()
labels_0 = [l.split()[0].strip('"') for l in new_data]
print('Labels:', Counter(labels_0))

for p, q in [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]:
    fig, ax = plt.subplots()
    ax.scatter(arr.T[p], arr.T[q])
    for i, txt in enumerate(labels):
        ax.annotate(txt, (arr.T[p, i], arr.T[q, i]))

dd = defaultdict(list)
dd = {}
ds = np.zeros((npop, npop))
stds = np.zeros((npop, npop))

ulabs = set(labs)
small = lambdas[:-1]
T = np.sqrt(2/np.average(small))

for i in ulabs:
    for j in ulabs:
        i_s = np.where(labs == i)
        j_s = np.where(labs == j)
        block = delta[i_s].T[j_s]*T/2
        block = block.reshape(block.shape[0]*block.shape[1],)
        if i == j:
            block = [el for el in block if el > 1e-5]
        ds[i, j] = np.average(block)
        stds[i, j] = np.std(block)
        dd[i, j] = block

plt.show()