import numpy as np
import pickle


def calc_mds(file: str, outfile: str) -> None:
    """Read distance matrix from 'file', calculate the eigensystem,
    and store it into 'outfile'.
    """
    with open(file, 'rb') as f:
         delta = pickle.load(f)

    N = len(delta)

    a = -delta**2/2

    at0 = np.sum(a, 0)/N
    att = np.sum(a)/N**2

    one = np.ones((N,))
    b = a - np.outer(at0, one) - np.outer(one, at0) + att

    lambdas, vecs = np.linalg.eigh(b)

    with open(outfile, 'wb') as f:
        pickle.dump((lambdas, vecs), f)