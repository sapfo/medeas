import numpy as np
import pickle

N = 294

with open('temp_asd.asd', 'rb') as f:
    delta = pickle.load(f)

a = -delta**2/2

at0 = np.sum(a, 0)/N
att = np.sum(a)/N**2

one = np.ones((N,))
b = a - np.outer(at0, one) - np.outer(one, at0) + att

lambdas, vecs = np.linalg.eig(b)

with open('temp_eig.data', 'wb') as f:
    pickle.dump((lambdas, vecs), f)