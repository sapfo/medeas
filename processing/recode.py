#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 14:28:40 2017

@author: ivan
"""

import numpy as np

def recode_wide(file: str, outfile: str) -> None:
    """Recode file in the "wide" format, for example:
        11101 -> 2 2 2 1 2
    """
    #with open(file) as f:
    #    data_lines = f.readlines()
    
    data = np.genfromtxt(file, dtype='int8', delimiter=1)
    np.savetxt(outfile, data + 1, fmt='%1d')
    print('Saved ', outfile)