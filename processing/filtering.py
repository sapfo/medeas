#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 16:27:49 2017

@author: ivan
"""

import numpy as np
from typing import List
import pickle
from options import FST




def set_missing(ancestry_file: str, infile: str, outfile: str,
                all_labels: List[str], groups: List[int]) -> None:
    """Set all data with ancestry from 'groups' as missing."""
    print("Processing file", ancestry_file)
    all_labels = np.array(all_labels)
    columns = np.where(all_labels == 'ABO')
    with open(infile, 'rb') as f:
        full_data = pickle.load(f)
    with open(ancestry_file, 'rb')as f:
        ancestry = pickle.load(f)
    if ancestry.shape[0]: # guard for empty files
        for group in groups:
            is_not_ancestry = np.sign(ancestry - group) ** 2
            full_data.T[columns] = full_data.T[columns] * is_not_ancestry.T
    with open(outfile, 'wb') as f:
        pickle.dump(full_data, f)


def filter_sparse(infile_pattern: str, outfile_pattern: str, ratio: float,
                  labels: List[str],chromosomes: range) -> 'np.ndarray[str]':
    """Keep only individuals where amount of non=missing sites is
    at least 'ratio'.
    """
    labels = np.array(labels)
    total = 0
    non_missing = None
    for n in chromosomes:
        print("First read: chromosome", n)
        with open(infile_pattern.format(n), 'rb') as f:
            data = pickle.load(f)
        total += data.shape[0]
        non_missing = (non_missing if non_missing is not None else
                       np.zeros((data.shape[1],))) + np.sum(np.sign(data), axis=0)
    columns = np.where(non_missing > ratio * total)
    for n in chromosomes:
        with open(infile_pattern.format(n), 'rb') as f:
            data = pickle.load(f)
        data = data.T[columns].T
        print("Writing file: chromosome", n)
        with open(outfile_pattern.format(n), 'wb') as f:
            pickle.dump(data, f)
    labels = labels[columns]
    return labels


def filter_manual(infile: str, outfile: str, pops: List[str],
                  labels: List[str]) -> 'np.ndarray[str]':
    """Manualy remove individuals from populations with labels 'pops'.

    'labels' is the initial list of labels for data. Returns the filtered list
    of labels the filtered SNP data is pickled.
    """
    short_labs = np.array([l.split()[0] for l in labels])
    print('Manual drop: ', pops, infile)
    labels = np.array(labels)
    with open(infile, 'rb') as f:
        data = pickle.load(f)

 #   if FST:
 #       compute_fst(short_labs,data)


    mask = np.ones((data.shape[1],))
    for pop in pops:
        mask[np.where(short_labs == pop)] = 0

    columns = np.where(mask)

    data = data.T[columns].T
    labels = labels[columns]
    with open(outfile, 'wb') as f:
        pickle.dump(data, f)
    return labels



def compute_fst(short_labs: 'np.array[str]',data: 'np.ndarray[int]'):
    where_BRI = np.where(short_labs == 'BRI')
    where_CHI = np.where(short_labs == 'CHI')
    data_BRI = data.T[where_BRI] - 1
    data_CHI = data.T[where_CHI] - 1

    len_BRI = len(short_labs[where_BRI])
    len_CHI = len(short_labs[where_CHI])
    print(len_BRI, len_CHI)

    p_BRI = np.mean(data_BRI, axis=0, dtype=np.float64)
    p_CHI = np.mean(data_CHI, axis=0, dtype=np.float64)
    # We take only variant sites
    non_zero = np.where((p_BRI != 0) | (p_CHI != 0))
    p_BRI = p_BRI[non_zero]
    p_CHI = p_CHI[non_zero]
    non_one = np.where((p_BRI != 1) | (p_CHI != 1))
    p_BRI = p_BRI[non_one]
    p_CHI = p_CHI[non_one]

    p = (p_BRI * len_BRI + p_CHI * len_CHI) / (len_BRI + len_CHI)
    combined = (p_BRI * (1 - p_BRI) * len_BRI +
                p_CHI * (1 - p_CHI) * len_CHI) / (len_BRI + len_CHI)
    total = p * (1 - p)
    print(f'F_ST = {np.mean(1 - combined/total)}')
    print(f'dF_ST = {np.std(1 - combined/total)}')