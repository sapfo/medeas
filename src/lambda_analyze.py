import numpy as np
from scipy.optimize import curve_fit
from typing import Sequence, Tuple
import matplotlib.pyplot as plt
import pickle

TESTING = True

def dens(x: float, a: float, b: float) -> float:
    """Integral of Marchenko-Pastur distribution function on interval
    a < x < b.
    """
    arc1 = np.arcsin((2*x - a - b)/(b - a))
    arc2 = np.arcsin(((a+b)*x-2*a*b)/x/(b - a))
    res = (np.sqrt((b-x)*(x-a))
            + (a+b)*arc1/2
            - np.sqrt(a*b)*arc2
            + np.pi*((a+b)/2 - np.sqrt(a*b))/2)
    return res


def find_TW(lambdas: Sequence[float], L: float) -> float:
    """Find Tracy-Widom statistics for the highest (first) eigenvalue
    among 'lambdas'.
    """
    m = len(lambdas)
    n = int(L)
    mu = (np.sqrt(n - 1) + np.sqrt(m))**2 / n
    sigma = (1/np.sqrt(n-1) + 1/np.sqrt(m))**(1/3) * (np.sqrt(n - 1) + np.sqrt(m)) / n
    l = m*lambdas[0]/np.sum(lambdas)
    print(l, mu, sigma)
    s = (l - mu) / sigma
    return s

def find_T_and_L(file: str) -> Tuple[float, float]:

    with open(file, 'rb') as f:
        lambdas, vecs = pickle.load(f)

    # finding L and T from fit
    lambdas_s = np.array(sorted(lambdas))
    N = len(lambdas)

    def dens_fit(x: float, T: float, L: int):
        """Integral of Marchenko-Pastur distribution function for
        arbitrary x from total tree length 'T' and number of markers L.
        """
        a = 2/T**2 * (1 - np.sqrt(N/L))**2
        b = 2/T**2 * (1 + np.sqrt(N/L))**2
        if x <= a + 0.0001:
            return 0.0
        if x >= b - 0.0001:
            return N * 1.0
        return N * dens(x, a, b) / dens(b-0.0001, a, b)

    # TODO: Fix this. Lambas should never be negative
    lambdas_s = np.array([max(l, 0.01) for l in lambdas_s])
    #lambdas_s = lambdas_s - lambdas_s.min() + 0.0001

    popt, pcov = curve_fit(np.vectorize(dens_fit, otypes=[np.float]),
                           lambdas_s, range(len(lambdas_s)),
                           p0 = (1, N+1), bounds=([0.1, N], [10, 100*N]))
    if TESTING:
        plt.plot(lambdas_s, range(len(lambdas_s)))
        lambdas_se = np.linspace(lambdas.min(), lambdas.max(), 5000)
        l_dens_fit = [dens_fit(l, *popt) for l in lambdas_se]
        plt.plot(lambdas_se, l_dens_fit)
        plt.show()
    print('FIT: ', popt, pcov)
    return popt

# ====================
# Tracy-Widom
# ====================

def find_K(file, L):

    with open(file, 'rb') as f:
        lambdas, vecs = pickle.load(f)

    lambdas0 = lambdas.copy()
    N = len(lambdas)

    s1 = np.sum(lambdas0)
    s2 = np.sum(lambdas0**2)

    print((N+1)*s1**2/((N-1)*s2 - s1**2))
    print(s1)
    print(s2)

    print('-'*20)
    s = 5
    i = 1
    while s > 3.2724:  # for p-value 0.001
        s = find_TW(lambdas, L)
        print(f'eigenvalue {i}: TW = {s}')
        i += 1
        lambdas = lambdas[1:]

    if TESTING:
        lambdas = list(sorted(lambdas, reverse=True))
        lambdas = np.array(lambdas)
        plt.figure()
        plt.plot(lambdas, 'b.')
        plt.figure()
        plt.hist(lambdas, 50)
    return i - 1
