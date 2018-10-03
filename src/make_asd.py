# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 15:29:37 2016

@author: ivan
"""

import numpy as np
import pickle


def  asd_main(pp: int, name: str, out_name: str, chromosomes: range,
             bootsize: int,
             txt_format: bool = False
             ) -> None:
    """Calculate the distance matrix with Minkowski parameter 'pp'.

    'name' is the input file name (or format string) with SNP data.
    'out_name' is the output binary file (pickle). If 'txt_format' is
    True, the SNP data will be read from a single text file (useful for
    tests with scrm), otherwise the data will be read from (binary)
    chromosome files.
    """

    with open(name) as f:
        print('Chunk loading started')
        N = len(f.readline()) // 2
        while True:
            print(f'N = {N}')
            data_lines = f.readlines()

            if not data_lines:
                break
            data = np.array([np.fromstring(line, sep=' ',  # line[-cut:-1]
                                           dtype='int8')
                             for line in data_lines])
        print('Chunk loading done')
    nb_loci = data.shape[0]
    for i in range(nb_loci):
        p = np.sum(data[i,:]-1)/N
        data[i,:] = data[i,:]/np.sqrt(p*(1-p))
    def compute_distance(data):
        delta = np.zeros((N,N))
        for i in range(N):
            for j in range(i + 1, N):
                delta[i, j] = delta[j, i] = sum(np.abs((data[:,i]-data[:,j])))**(1/pp)/N
        return(delta)

    delta = compute_distance(data)
    with open(out_name, 'wb') as f:
        pickle.dump(delta, f)


