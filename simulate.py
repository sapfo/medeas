#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 17:03:27 2017

@author: ivan
"""

from subprocess import Popen, PIPE
import sys

def run_simulation_two_pops(n1: int, n2: int, L: int, theta: float, D: float):

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
    
run_simulation_two_pops(30, 70, 100, 1, 0.1)