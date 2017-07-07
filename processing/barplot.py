#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 19:30:31 2017

@author: ivan
"""

import numpy as np
import matplotlib.pyplot as plt

total_bars = {}
total_freqs = {}

zero = np.zeros((118,))
for group in range(1, 5):
    total_bars[group] = np.loadtxt(f'abo.all.bars.{group}.txt')

print('Plotting........')

bins = np.arange(0.5, 119.5, 1)
colors = ['purple', 'red', 'yellow', 'blue']
colors = ['yellow', 'blue', 'red', 'purple']

# labels = np.loadtxt

plt.figure()
width = 1
for i, group in enumerate([4, 3, 1, 2]):
    plt.bar(np.arange(118), total_bars[group], width=width,
            bottom = sum((total_bars[g] for g in [4, 3, 1, 2][:i]), zero),
            color=colors[group-1])

plt.show()