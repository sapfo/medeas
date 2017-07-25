import numpy as np
import pickle

#N = 294

import sys

def calc_mds(file, outfile):

    with open(file, 'rb') as f: # temp_asd.asd
         delta = pickle.load(f)

    N = len(delta)
    print(N)

    a = -delta**2/2

    at0 = np.sum(a, 0)/N
    att = np.sum(a)/N**2

    one = np.ones((N,))
    b = a - np.outer(at0, one) - np.outer(one, at0) + att

    lambdas, vecs = np.linalg.eig(b)

    with open(outfile, 'wb') as f:
        pickle.dump((lambdas, vecs), f)