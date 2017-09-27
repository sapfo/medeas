#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 17:03:27 2017

@author: ivan
"""

from subprocess import Popen, PIPE
import sys
import numpy as np

TWO_POPS = True

def run_simulation_two_pops(n1: int, n2: int, L: int, theta: float, D: float):

    with open('res8.txt', 'a') as f:
        f.write(f'{L} {2*D}')

    scrm = Popen(f'../scrm/scrm {n1+n2} {L} -t {theta} -I 2 {n1} {n2} -ej {D} 1 2'
                 ' --print-model -l -1 -L'.split(' '), stdout=PIPE)
    grep = Popen(['grep',  rb'^\d*$'], stdin=scrm.stdout, stdout=PIPE)
    scrm.stdout.close()
    output = grep.communicate()
    data = output[0].decode('utf-8')

    with open('../scrm/tmpout.txt', 'w') as f:
        f.write(data)

    trans = Popen('python transcode.py ../scrm/tmpout.txt ../scrm/tmpres.txt'.split(' '),
                  stdout=sys.stdout)
    trans.communicate()

    medeas = Popen(['time',
                    'python',
                    'main.py',
                    '../raw_data_copy/eur.chi.pap.wcd.abo.chr{}.g10.txt.0.Viterbi.txt',
                    '../raw_data_copy/eur.chi.pap.wcd.abo.chr{}.txt',
                    '../raw_data_copy/haplotype_labels.txt',
                    '-asd',
                    '-analyze'], stdout=sys.stdout)
    medeas.communicate()

if TWO_POPS:

    Ls = np.logspace(2, 5, num=100)
    Ls2 = [300, 1000]

    for D in [0.02, 0.01]:
        for L in Ls2:
            for _ in range(10):
                run_simulation_two_pops(30, 70, int(L), 1, D)

def run_simulation_three_pops(n1: int, n2: int, n3: int, L: int, theta: float, D: float, D1: float):

    with open('res3.txt', 'a') as f:
        f.write(str(L))
    
    scrm = Popen(f'../scrm/scrm {n1+n2+n3} {L} -t {theta} -I 3 {n1} {n2} {n3} -ej {D} 2 1 -ej {D1} 3 2'
                 ' --print-model -l -1 -L'.split(' '), stdout=PIPE)
    grep = Popen(['grep',  rb'^\d*$'], stdin=scrm.stdout, stdout=PIPE)
    scrm.stdout.close()
    output = grep.communicate()
    data = output[0].decode('utf-8')

    with open('../scrm/tmpout.txt', 'w') as f:
        f.write(data)

    trans = Popen('python transcode.py ../scrm/tmpout.txt ../scrm/tmpres.txt'.split(' '),
                  stdout=sys.stdout)
    trans.communicate()

    medeas = Popen(['time',
                    'python',
                    'main.py',
                    '../raw_data_copy/eur.chi.pap.wcd.abo.chr{}.g10.txt.0.Viterbi.txt',
                    '../raw_data_copy/eur.chi.pap.wcd.abo.chr{}.txt',
                    '../raw_data_copy/haplotype_labels.txt',
                    '-asd',
                    '-analyze'], stdout=sys.stdout)
    medeas.communicate()

if not TWO_POPS:

    Ls3 = [100, 300, 1000, 3000, 10000, 30000]
    for L in Ls3:
        for _ in range(10):
            run_simulation_three_pops(20, 30, 50, int(L), 1, 0.2, 0.1)