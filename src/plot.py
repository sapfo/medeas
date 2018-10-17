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
import matplotlib.lines as mlines
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
    markers = [".","v","*","+","x","2","p","^","s"]
    for p in range(len(np.unique(labels_inferred))+2):
        for q in range(p+1,len(np.unique(labels_inferred))+2):
            fig, ax = plt.subplots(figsize=(15,10))
            for population_index, population_name in enumerate(np.unique(label_given)):
                position_population = np.where(population_name == label_given)
                index_colors = labels_inferred[position_population]
                color_value = [colors[index_color] for index_color in index_colors]
                markers_value = markers[population_index]
                ax.scatter(arr.T[p,position_population].ravel(), arr.T[q,position_population].ravel(), c = color_value ,marker = markers_value,s=100)
            plt.legend(np.unique(label_given))
            leg = ax.get_legend()
            for point in leg.legendHandles:
                point.set_color('black')
            dir_plot = os.path.join(output_folder,"mds_plot")
            if not os.path.isdir(dir_plot):
                os.mkdir(dir_plot)
            markers_shape = [mlines.Line2D([], [], color='black', marker=marker_shape, linestyle='None') for marker_shape in markers]
            markers_color = [mlines.Line2D([], [], color=marker_color, marker="o", linestyle='None') for marker_color in colors]
            legend1 = plt.legend(markers_shape,np.unique(label_given) , loc=1,title="True population")
            legend2 = plt.legend(markers_color,np.unique(labels_inferred), loc=2,title="Inferred population")
            ax.add_artist(legend1)
            ax.add_artist(legend2)
            ax.set_xlabel(f'PC. {p+1}')
            ax.set_ylabel(f'PC. {q+1}')
            fig.savefig(os.path.join(dir_plot, f'mds_axis{p+1}_{q+1}.pdf'))


