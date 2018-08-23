#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 17:03:27 2017

@author: ivan
"""
import yappi
from subprocess import Popen, PIPE
import sys
import os
import datetime

yappi.start()

folder = sys.argv[1] # folder should exist. its currently created by job_sim.sh

TWO_POPS = False
FNAME = 'res_sim.txt'
scrm_file_raw_result = os.path.join(folder, '_tmp_out.txt')
scrm_file_clean_result = os.path.join(folder, '_tmp_res.txt')
scrm_file_seed = os.path.join(folder, '_tmp_seed.txt')
file_fake_labs = os.path.join(folder, 'fake_labs.txt')


def run_simulation_two_pops(n1: int, n2: int, L: int, theta: float, D: float):
    scrm_exec = sys.argv[2]
    # n1 is always an outgroup

    # create fake labels because medeas needs that
    with open(file_fake_labs, 'w') as f:
        for i in range(n1+n2):
            if i < n1:
                pop = 'PAP'
            else:
                pop = 'BRI'
            f.write(f'{pop}  {pop}{i}\n')

    # run scrm  for simulations
    scrm_command = f'{scrm_exec} {n1+n2} {L} -t {theta} -I 2 {n1} {n2} -ej {D} 1 2 --print-model -l -1 -L'
    with open(os.path.join(folder, 'scrm_command.txt'), 'w') as f:
        f.write(scrm_command)
    scrm = Popen(scrm_command.split(' '), stdout=PIPE)
    grep = Popen(['grep',  rb'^[0-9]*$'], stdin=scrm.stdout, stdout=PIPE)
    scrm.stdout.close()
    output = grep.communicate()
    data = output[0].decode('utf-8')

    # get the data from scrm, incl seed

    with open(scrm_file_raw_result, 'w') as f:
        f.write(data)

    trans = Popen(f'python transcode.py {scrm_file_raw_result} {scrm_file_clean_result} {scrm_file_seed} {n1+n2}'.split(' '),
                  stdout=sys.stdout)
    trans.communicate()

    #read and write the seed
    with open(scrm_file_seed, 'r') as f:
        seed = f.read().strip()
    with open(os.path.join(folder, FNAME), 'a') as f:
        f.write(f'{seed} {L} {2*D}') #RESCALE THE SPLIT TIME

    #run medeas, that is main.py, on the simulated data
    # simulation flag is hard coded for now: SIMULATION = True, needs changing
    #arguments for medeas:
    #argv[1] : ancestry pattern (painting) -- never read, SIMULATION = True
    #argv[2] : SNP data -- never read, SIMULATION = True
    #argv[3] : labels file -- neer read, SIMULATION = True
    #argv[4:-4]: flags, -asd: calculate the matrix, -analyze: run all analyses
    #argv[-4]: true nssites
    #argv[-3]: true split time
    #argv[-2]: folder where the simulated data lives
    #argv[-1]: filename where to store the summary of the simulation
    medeas = Popen(['python','main.py',
                    'TEMP',
                    'TEMP',
                    'TEMP',
                    '-asd', '-analyze', str(L), str(2*D), folder, FNAME])
    medeas.communicate()

if TWO_POPS:
    # run for a range lf nssites
    Ls = [10**i for i in range(2, 5)]
    with open(os.path.join(folder, FNAME), 'a') as f:
        f.write(f'Seed L D D_inf D_std K_inf T L_inf\n')

    # medeas is run on simulation mode so it will append the inferred values
    D = float(sys.argv[3])
    for L in Ls:
        run_simulation_two_pops(50, 150, int(L), 1, D)

def run_simulation_three_pops(n1: int, n2: int, n3: int, L: int, theta: float, D: float, D1: float):
    # DOES NOT WORK YET
    # n1 is always an outgroup

    with open(file_fake_labs, 'w') as f:
        for i in range(n1+n2+n3):
            if i < n1:
                pop = 'PAP'
            elif i < n2:
                pop = 'BRI'
            else:
                pop = 'CHI'
            f.write(f'{pop}  {pop}{i}\n')

    with open('res3.txt', 'a') as f:
        f.write(str(L))
    
    scrm = Popen(f'scrm {n1+n2+n3} {L} -t {theta} -I 3 {n1} {n2} {n3} -ej {D} 2 1 -ej {D1} 3 2'
                 ' --print-model -l -1 -L'.split(' '), stdout=PIPE)
    grep = Popen(['grep',  rb'^[0-9]*$'], stdin=scrm.stdout, stdout=PIPE)
    scrm.stdout.close()
    output = grep.communicate()
    data = output[0].decode('utf-8')

    with open(scrm_file_raw_result, 'w') as f:
        f.write(data)

    trans = Popen(f'python transcode.py {scrm_file_raw_result} {scrm_file_clean_result} {scrm_file_seed} {n1+n2+n3}'.split(' '),
                  stdout=sys.stdout)
    trans.communicate()
    medeas = Popen(['python','main.py',
                    'TEMP',
                    'TEMP',
                    'TEMP',
                    '-asd', '-analyze', str(L), str(2*D), folder, FNAME])
    medeas.communicate()

if not TWO_POPS:
    Ls3 = [10000]
    for L in Ls3:
        run_simulation_three_pops(50, 100, 150, int(L), 1, 0.2, 0.1)




func_stats = yappi.get_func_stats()
#func_stats.save('callgrind.out.' + datetime.datetime.now().isoformat(), 'CALLGRIND')
yappi.stop()
yappi.clear_stats()