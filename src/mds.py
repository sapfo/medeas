import numpy as np


a = -delta**2/2

at0 = np.sum(a, 0)/N
att = np.sum(a)/N**2

one = np.ones((N,))
b = a - np.outer(at0, one) - np.outer(one, at0) + att

lambdas, vecs = np.linalg.eig(b)