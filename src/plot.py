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

import os

def plot_mds(arr, label_given, labels_inferred,output_folder):
    """Plot the MDS plot
    """
    # TODO: autogenerate nice summary plot depending on 'npop'
    label_given = np.array(label_given)
    label_given_index = np.copy(label_given)
    for index_label, label in enumerate(np.unique(label_given)):
        label_given_index[label_given == label] = index_label
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    markers = [".","s","v","*","+","x","2","p","^"]
    #for p, q in [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]:
    for p in range(len(np.unique(labels_inferred))):
        for q in range(p+1,len(np.unique(labels_inferred))):
            fig, ax = plt.subplots()
            for population_index, population_name in enumerate(np.unique(label_given)):
                position_population = np.where(population_name == label_given)
                index_colors = labels_inferred[position_population]
                color_value = [colors[index_color] for index_color in index_colors]
                markers_value = markers[population_index]
                ax.scatter(arr.T[p,position_population].ravel(), arr.T[q,position_population].ravel(), c = color_value ,marker = markers_value,s=100)
            fig.savefig(os.path.join(output_folder, f'mds_axis{p+1}_{q+1}.pdf'))


