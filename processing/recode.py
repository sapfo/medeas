#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 14:28:40 2017

@author: ivan
"""

import numpy as np
import pickle

def recode_wide(file: str, outfile: str,
                anc_file: str, out_anc_file: str) -> None:
    """Recode file in the "wide" format, for example:
        11101 -> 2 2 2 1 2
    """
    #with open(file) as f:
    #    data_lines = f.readlines()
    
    data = np.genfromtxt(file, dtype='int8', delimiter=1)
    with open(outfile, 'wb') as f:
        pickle.dump(data + 1, f)
    anc_data = np.genfromtxt(anc_file, dtype='int8')
    with open(out_anc_file, 'wb') as f:
        pickle.dump(anc_data, f)
    print('Saved ', out_anc_file)