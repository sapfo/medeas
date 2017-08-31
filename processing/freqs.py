#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 11:31:18 2017

@author: ivan
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

from typing import Tuple
from collections import defaultdict
import pickle

freqs = defaultdict(list)
bars = defaultdict(list)

def process_chromosome(file: str) -> None:
    print(f'Started processing', file)
    
    with open(file, 'rb') as f:
        ancestry = pickle.load(f)
    print(f'Loaded ancestry from', file)

    group: int
    for group in range(1, 5):
        is_ancestry = 1 - np.sign(ancestry - group) ** 2
        print(f'Processed selectors for group {group}')
        freqs[group].append(np.sum(is_ancestry, axis=1))
        print(f'Processed frequencies for group {group}')
        bars[group].append(np.sum(is_ancestry, axis=0))
        print(f'Processed bars for group {group}')

def calculate_freqs(pattern, outpattern):
    for num in range(1, 23):
        process_chromosome(pattern.format(num))

    total_bars = {}
    total_freqs = {}

    zero = np.zeros((118,))
    for group in range(1, 5):
        total_bars[group] = sum(bars[group], zero)
        total_freqs[group] = np.concatenate(tuple(freqs[group]))
        with open(outpattern.format(group) + '.bars', 'wb') as f:
            pickle.dump(total_bars[group], f)
        with open(outpattern.format(group) + '.freqs', 'wb') as f:
            pickle.dump(total_freqs[group], f)


if __name__ == '__main__':
    print('Plotting........')

    bins = np.arange(0.5, 119.5, 1)
    colors = ['black', 'blue', 'green', 'red']
    for group in range(1, 5):
        plt.hist(total_freqs[group], bins=bins, color=colors[group-1])

    plt.figure()
    width = 0.5
    for group in range(1, 5):
        plt.bar(np.arange(118), total_bars[group], width=width,
                bottom = sum((total_bars[g] for g in range(1, group)), zero),
                color=colors[group-1])

    plt.show()