# -*- coding: utf-8 -*-
"""
Created on Thu Sep 20 15:29:37 2018

@author: Frederic
"""




import os

from processing.recode import recode_wide
from processing.freqs import calculate_freqs
from processing.filter_freqs import load_freqs, soft_filter
from processing.hardfilter import hard_filter
from processing.filtering import set_missing, filter_sparse, filter_manual
from concurrent.futures import ProcessPoolExecutor as PPE
from multiprocessing import cpu_count


def preprocess_data(snps_pattern: str, ancestry_pattern: str, output_folder:str,output_file: str, chromosomes:range):
    snps_pattern_stped = os.path.join(snps_pattern + '.stped')
    group_pattern = os.path.join(output_folder, 'group.{}')
    freqs_pattern = group_pattern + '.freqs'
    filtered_pattern = snps_pattern_stped + '.filtered'
    filtered_ancestry_pattern = ancestry_pattern + '.filtered'
    hard_filtered_pattern = filtered_pattern + '.hard'
    hard_filtered_ancestry_pattern = filtered_ancestry_pattern + '.hard'

    missingnes_pattern = os.path.join(output_folder, 'processed', 'chr.{}.stped')
    directory = os.path.join(output_folder, 'processed')
    if not os.path.exists(directory):
        os.makedirs(directory)


    def do(task: Callable[..., None], count: int = cpu_count()):
        """A simple helper to paralelize given task across chromosomes."""
        executor = PPE(count)
        return list(executor.map(task, chromosomes))


    def read_labs(file: str) -> List[str]:
        with open(file) as f:
            return f.readlines()


    def save_labs(labs: Iterable[str], file: str) -> None:
        with open(file, 'w') as f:
            return f.writelines(labs)

    def recode(n):
        recode_wide(snps_pattern.format(n), snps_pattern_stped.format(n),
                    ancestry_pattern.format(n),
                    ancestry_pattern.format(n) + '.recoded')
    do(recode, 6)


    calculate_freqs(ancestry_pattern + '.recoded', group_pattern, chromosomes)


    for n in chromosomes:
        new_labs = filter_manual(snps_pattern_stped.format(n),
                                 snps_pattern_stped.format(n) + '.selected',
                                 ["CHI", "BRI"], read_labs(labels_file))
    save_labs(new_labs, labels_file + '.selected')


    mu, sigma = load_freqs(freqs_pattern)


    def filterf(n):
        soft_filter(ancestry_pattern.format(n) + '.recoded',
                    snps_pattern_stped.format(n) + '.selected',
                    filtered_pattern.format(n),
                    filtered_ancestry_pattern.format(n), mu, sigma, 5),

    do(filterf)


    def hardfilt(n):
        hard_filter(filtered_ancestry_pattern.format(n),
                    hard_filtered_ancestry_pattern.format(n),
                    filtered_pattern.format(n),
                    hard_filtered_pattern.format(n),
                    [4, 3, 2, 1], 0.1)  # PAP or WCD

    do(hardfilt)


    labs = np.array([l.split()[0] for l in read_labs(labels_file)])

    def setmiss(n):
        set_missing(hard_filtered_ancestry_pattern.format(n),
                    hard_filtered_pattern.format(n),
                    output_file.format(n),
                    labs, [])


    do(setmiss)

    new_labs = filter_sparse(missingnes_pattern, missingnes_pattern + '.filtered',
                             0.3, read_labs(labels_file + '.selected'), chromosomes)
    save_labs(new_labs, labels_file + '.filtered')