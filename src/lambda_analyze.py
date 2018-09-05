import numpy as np
from scipy.optimize import curve_fit
from typing import Sequence, Tuple
import matplotlib.pyplot as plt
import pickle

from options import TESTING


def dens(x: float, a: float, b: float) -> float:
    """Integral of Marchenko-Pastur distribution function on interval
    a < x < b.
    """
    arc1 = np.arcsin((2*x - a - b)/(b - a)/1.00001)
    arc2 = np.arcsin(((a+b)*x-2*a*b)/x/(b - a)/1.00001)
    res = (np.sqrt((b-x)*(x-a))
            + (a+b)*arc1/2
            - np.sqrt(a*b)*arc2
            + np.pi*((a+b)/2 - np.sqrt(a*b))/2)
    return res


def find_TW(lambdas: Sequence[float], L: float, T: float, N: int) -> float:
    """Find Tracy-Widom statistics for the highest (first) eigenvalue
    among 'lambdas'.
    """
    m = len(lambdas)
    n = int(L)
    mu = (np.sqrt(n - 1) + np.sqrt(m))**2 / n
    sigma = (1/np.sqrt(n-1) + 1/np.sqrt(m))**(1/3) * (np.sqrt(n - 1) + np.sqrt(m)) / n
    l = lambdas[0] * T**2/2 # * m/np.sum(lambdas)
    s = (l - mu) / sigma
    print(f'n = {n}, l = {l},mu =  {mu},sigma = {sigma},T**2/2 = {T**2/2}, m/np.sum(lambdas) = {m/np.sum(lambdas)} s = {s}')
    return s


def find_T_and_L(file: str) -> Tuple[float, float]:
    """Find total tree length T and effective number of markers L using
    the bulk eigenvalues from eigensystem stored in (pickled) 'file'.
    """
    with open(file, 'rb') as f:
        lambdas, vecs = pickle.load(f)

    # finding L and T from fit

    lambdas_s = np.array(sorted(lambdas))
    lambdas_s = lambdas_s[1:]
    N = len(lambdas)

    def dens_fit(x: float, T: float, L: int):
        """Integral of Marchenko-Pastur distribution function for
        arbitrary x from total tree length 'T' and number of markers L.
        """
        a = 2/T**2 * (1 - np.sqrt(N/L))**2
        b = 2/T**2 * (1 + np.sqrt(N/L))**2
        if x <= a:
            return 0.0
        if x >= b:
            return (N-1) * 1.0
        return (N-1) * dens(x, a, b) / dens(b, a, b)

    popt, pcov = curve_fit(np.vectorize(dens_fit, otypes=[np.float]),
                           lambdas_s, range(len(lambdas_s)),
                           p0 = (1, N+1), bounds=([0.1, N], [100, 1000*N]))

    # popt[1] = 2000
    if TESTING:
        plt.close()
        plt.plot(lambdas_s, range(len(lambdas_s)))
        lambdas_se = np.linspace(lambdas.min(), lambdas.max(), 5000)
        l_dens_fit = [dens_fit(l, *popt) for l in lambdas_se]
        import matplotlib.pylab as pyl
        pyl.xlim([0.0,0.2])
        plt.plot(lambdas_se, l_dens_fit)
        plt.savefig('fit.pdf')
        plt.clf()
        plt.hist(lambdas_s)
        plt.savefig('fit2.pdf')
        plt.clf()
        plt.hist(l_dens_fit)
        plt.savefig('fit3.pdf')
    p_err = np.sqrt(np.diag(pcov))
    print('FIT: ', popt, p_err)
    return popt[0], popt[1] - p_err[1]  # Take L 1σ lower to compensate
                                        # for tendency to overestimate

# ====================
# Tracy-Widom
# ====================

def find_K(file: str, L: float, T: float) -> int:
    """Find the number of clusters K using the Tracy-Widom statistics for
    large eigenvalues. 'L' is the effective number of markers, 'T' is
    the total tree length, the eigensystem is read from file named 'file'.
    """
    with open(file, 'rb') as f:
        lambdas, vecs = pickle.load(f)

    N = len(lambdas)
    lambdas = list(sorted(lambdas, reverse=True))
    lambdas0 = np.array(lambdas)

    s1 = np.sum(lambdas0)
    s2 = np.sum(lambdas0**2)
    L_patterson =  (N+1)*s1**2/((N-1)*s2 - s1**2)

    print("Used estimator for L ", L)
    print("Patterson estimator for L (not used)",L_patterson)  # Patterson estimator for L

    print('-'*20)
    s = 5
    i = 1
    while s > 3.2724:  # for p-value 0.001
        s1 = np.sum(np.array(lambdas))
        s2 = np.sum(np.array(lambdas) ** 2)
        L_patterson = (N + 1) * s1 ** 2 / ((N - 1) * s2 - s1 ** 2)
        print("Patterson estimator for L (not used)", L_patterson)  # Patterson estimator for L

        s = find_TW(lambdas, L, T, N)
        print(f'eigenvalue {i}: TW = {s}')
        i += 1
        lambdas = lambdas[1:]

    if TESTING:
        lambdas = list(sorted(lambdas0, reverse=True))
        lambdas = np.array(lambdas)
        plt.figure()
        plt.plot(lambdas, 'b.')
        plt.figure()
        plt.hist(lambdas, 50)
    return i - 1
