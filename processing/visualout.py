#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:10:15 2017

@author: ivan
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

from collections import defaultdict

wcd_ancestry = np.loadtxt('abo.all.bars.4.txt')

chunk_size = 1000
samples = np.where(wcd_ancestry < 2e6)[0]

with open('haplotype_labels.txt') as f:
    labels0 = f.readlines()
labels = np.array([l.split()[1] for l in labels0])

groups = np.array([l.split()[0] for l in labels0])

samples = np.array([8, 9, 26, 27, 36, 37, 72, 73, 82, 83])

labels = labels[np.where(groups == 'ABO')]

print(labels)
print(groups)

print(labels[samples])

samples = np.concatenate((np.array([29]), samples))

def process_chromosome(num):

    ancestries = defaultdict(list)    
    ancestry = np.genfromtxt(f'eur.chi.pap.wcd.abo.chr{num}.g10.txt.0.Viterbi.txt',
                             dtype='int8')
    for index in samples:
        ancestries[index].append(ancestry.T[index])
        print(f'Processed chromosome #{num}, outlier #{index}')
    
    ancestry_total = {}
    ancestry_smooth = {}
    for index in samples:
        ancestry_total[index] = np.concatenate(tuple(ancestries[index]))
        ancestry_smooth[index] = np.zeros((ancestry_total[index].shape[0]//chunk_size,))
        for i in range(ancestry_smooth[index].shape[0]):
            ancestry_smooth[index][i:i+1] = np.sum(ancestry_total[index][chunk_size*i:
                                                                         chunk_size*i+chunk_size],
                                                   axis=0)
    
    aspect = 200
    plt.figure()
    f, axarr = plt.subplots(len(samples), sharex=True)
    for i, index in enumerate(samples):
        axarr[i].get_yaxis().set_ticks([])
        axarr[i].set_aspect(aspect)
        axarr[i].set_ylim([0, 1])
        axarr[i].pcolormesh(ancestry_smooth[index].reshape(1, ancestry_smooth[index].shape[0]),
                            vmin=chunk_size, vmax=4*chunk_size)
    
    axarr[len(samples)-1].set_xlabel(f'Site/{chunk_size}')    
    axarr[len(samples)//2].set_ylabel(f'haplotypes {list(reversed(samples))}')
    axarr[0].set_title(f'Chromosome #{num}')
    plt.savefig(f'Ancestry.outliers.{num}.pdf')
    print('Figure saved')

for num in range(int(sys.argv[1]), int(sys.argv[2])):
    process_chromosome(num)