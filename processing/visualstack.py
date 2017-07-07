#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 18:35:25 2017

@author: ivan
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

with open('haplotype_labels.txt') as f:
    labels0 = f.readlines()
labels = np.array([l.split()[1] for l in labels0])
groups = np.array([l.split()[0] for l in labels0])
labels = labels[np.where(groups == 'ABO')]

chunk_size = 1000
samples = int(sys.argv[1]), int(sys.argv[2])
chromosomes = range(1, 23)

sample_labels = labels[np.array(samples)]

aspect = 200
f, axarr = plt.subplots(2*len(chromosomes), sharex=True)

for order, num in enumerate(chromosomes):
    ancestry = np.genfromtxt(f'eur.chi.pap.wcd.abo.chr{num}.g10.txt.0.Viterbi.txt',
                             dtype='int8')
    for subord, index in enumerate(samples):

        print(f'Processed chromosome #{num}, outlier #{index}')
 
        ancestry_smooth = np.zeros((ancestry.T[index].shape[0]//chunk_size,))
        for i in range(ancestry_smooth.shape[0]):
            ancestry_smooth[i:i+1] = np.sum(ancestry.T[index][chunk_size*i:
                                                              chunk_size*i+chunk_size],
                                                   axis=0)
        i = 2*order + subord
        axarr[i].get_yaxis().set_ticks([])
        axarr[i].set_aspect(aspect)
        axarr[i].set_ylim([0, 1])
        axarr[i].pcolormesh(ancestry_smooth.reshape(1, ancestry_smooth.shape[0]),
                            vmin=chunk_size, vmax=4*chunk_size)

axarr[2*len(chromosomes)-1].set_xlabel(f'Site/{chunk_size}')    
axarr[len(chromosomes)].set_ylabel('Chromosomes')
axarr[0].set_title(f'Haplogroups {sample_labels}')
plt.savefig(f'Ancestry.outliers.samples.{samples[0]}.{samples[1]}.pdf')
print('Figure saved')
