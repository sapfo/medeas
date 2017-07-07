#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 13:30:43 2017

@author: ivan
"""

import numpy as np
import pickle

with open('austr.pp.2.asd', 'rb') as f:
    delta = pickle.load(f)

with open('../test/haplotype_labels.txt') as f:
    lines = f.readlines()

labels = np.array([l.split()[0] for l in lines])

print(labels)

where = np.where(np.logical_or(labels == 'ABO', labels == 'WCD'))[0]

print(where)

with open('austr.pp.2.cut.asd', 'wb') as f:
    pickle.dump(delta[where].T[where], f)