#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 17:03:27 2017

@author: ivan
"""

from subprocess import Popen, PIPE
import sys
import os
import datetime



folder = sys.argv[1] # folder should exist. its currently created by job_sim.sh

SIMULATION_CASE = 6
# 1 = 1 population
# 2 = 2 population
# 2 = 3 pop
# > 3 = n pop
FNAME = 'res_sim.txt'
scrm_file_raw_result = os.path.join(folder, '_tmp_out.txt')
scrm_file_clean_result = os.path.join(folder, '_tmp_res.txt')
scrm_file_seed = os.path.join(folder, '_tmp_seed.txt')
file_fake_labs = os.path.join(folder, 'fake_labs.txt')
scrm_exec = sys.argv[2]

def run_scrm(scrm_command: str):
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

def transcode_scrm(n: int):
    trans = Popen(f'python transcode.py {scrm_file_raw_result} {scrm_file_clean_result} {scrm_file_seed} {n}'.split(' '),
                  stdout=sys.stdout)
    trans.communicate()

    #read and write the seed
    with open(scrm_file_seed, 'r') as f:
        seed = f.read().strip()
    with open(os.path.join(folder, FNAME), 'a') as f:
        f.write(f'{seed}')

def run_simulation_single_pop(nb_individual: int, L: int, theta: float):
    """Run a simulation for nb_pop with each n individual, splitting at time D, 2D, 3D, ..."""
    with open(file_fake_labs, 'w') as f:
        for i in range(nb_individual):
            pop = 'pop'
            f.write(f'{pop} {pop}{i}\n')


    scrm_command = f'scrm {nb_individual} {L} -t {theta} --print-model -l -1 -L'
    run_scrm(scrm_command)
    transcode_scrm(nb_individual)

if SIMULATION_CASE == 1:
    Ls3 = [4000]
    for L in Ls3:
        run_simulation_single_pop(400, int(L), 1)

def run_simulation_two_pops(n1: int, n2: int, L: int, theta: float, D: float):
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
    run_scrm(scrm_command)
    transcode_scrm(n1+n2)

if SIMULATION_CASE == 2:
    Ls = [2000]
    D = 8
    for L in Ls:
        run_simulation_two_pops(100, 100, int(L), 2, D)

def run_simulation_three_pops(n1: int, n2: int, n3: int, L: int, theta: float, D: float, D1: float):
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
    scrm_command = f'scrm {n1+n2+n3} {L} -t {theta} -I 3 {n1} {n2} {n3} -ej {D} 2 1 -ej {D1} 3 2'
    run_scrm(scrm_command)
    transcode_scrm(n1+n2)



if SIMULATION_CASE == 3:
    Ls3 = [10000]
    for L in Ls3:
        run_simulation_three_pops(50, 100, 150, int(L), 1, 0.2, 0.1)


def run_simulation_n_pops(n_per_pop: int, nb_pop: int, L: int, theta: float, D: float):
    """Run a simulation for nb_pop with each n individual, splitting at time D, 2D, 3D, ..."""
    with open(file_fake_labs, 'w') as f:
        for i in range(n_per_pop * nb_pop):
            pop = 'pop{}'.format(i // n_per_pop)
            f.write(f'{pop} {pop}{i}\n')

    all_pop_string = nb_pop*(str(n_per_pop) + " ")
    all_split_string = " ".join([f'-ej {D*i} {i} {i+1}' for i in range(1,nb_pop)])
    scrm_command = f'scrm {n_per_pop*nb_pop} {L} -t {theta} -I {nb_pop} {all_pop_string}{all_split_string} --print-model -l -1 -L'
    run_scrm(scrm_command)
    transcode_scrm(n_per_pop * nb_pop)

if SIMULATION_CASE > 3:
    Ls3 = [10000]
    for L in Ls3:
        run_simulation_n_pops(50, SIMULATION_CASE, int(L), 1, 0.5)



def run_medeas():
    #run medeas, that is main.py, on the simulated data
    medeas = Popen(['python','main.py',
                    '-asd', '-analyze', folder, FNAME])
    medeas.communicate()

run_medeas()