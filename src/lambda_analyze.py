import numpy as np
from scipy.optimize import curve_fit
from typing import Sequence, Tuple
import matplotlib.pyplot as plt
import pickle





def find_T_and_t_within(file: str, label, populations) -> Tuple[float, float]:
    """Find total tree length T and effective number of markers L using
    the bulk eigenvalues from eigensystem stored in (pickled) 'file'.
    """
    with open(file, 'rb') as f:
        delta = pickle.load(f)

    np_population = len(populations)
    ts_over_T = np.zeros(np_population)
    for index_population, population in enumerate(populations):
        position_pop = np.where(label == population)[0]
        nb_individual = len(position_pop)
        sub_matrix = delta[np.ix_(position_pop,position_pop)]
        upper_sub_matrix = np.triu_indices(nb_individual, k = 1)
        t_over_T = 0.5*np.mean(sub_matrix[upper_sub_matrix])
        ts_over_T[index_population] = t_over_T
    print(f"ts/T =  {ts_over_T}")
    first_t = ts_over_T[0]
    T = 1/first_t
    ts = ts_over_T/first_t
    return((T,ts))

