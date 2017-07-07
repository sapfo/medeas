#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 18:03:48 2017

@author: ivan
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

num = int(sys.argv[1])
chunk_size = 100
samples = 118

ancestry = np.genfromtxt(f'eur.chi.pap.wcd.abo.chr{num}.g10.txt.0.Viterbi.txt',
                         dtype='int8')

ancestry_smooth = np.zeros((ancestry.shape[0]//chunk_size, ancestry.shape[1]))

for i in range(ancestry_smooth.shape[0]):
    ancestry_smooth[i:i+1] = np.sum(ancestry[chunk_size*i: chunk_size*i+chunk_size],
                                    axis=0)

aspect = 20

f, axarr = plt.subplots(samples, sharex=True)
for i in range(samples):
    axarr[i].get_yaxis().set_ticks([])
    axarr[i].set_aspect(aspect)
    axarr[i].set_ylim([0, 1])
    axarr[i].pcolormesh(ancestry_smooth.T[i:i+1], vmin=chunk_size, vmax=4*chunk_size)

axarr[samples-1].set_xlabel(f'Site/{chunk_size}')    
axarr[samples//2].set_ylabel(f'Individulas ({samples} haplotypes)')
axarr[0].set_title(f'Chromosome #{num}')
plt.savefig(f'Ancestry{num}.pdf')
print('Figure saved')

plt.figure()
plt.pcolormesh(ancestry_smooth, vmin=chunk_size, vmax=4*chunk_size)
plt.colorbar()

plt.show()
