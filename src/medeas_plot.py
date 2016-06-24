#!/usr/bin/env python

# Copyright (C) 2015,2016  Ole Tange, Mike DeGiorgio, Anna-Sapfo
# Malaspinas, Jose Victor Moreno-Mayar, Yong Wang, Ivan Levkivskyi and Free Software
# Foundation, Inc.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License,
# or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# For Interative use:
#    eval(parse(text=readLines("dist/src/bammds_plot.r")))
# sort( sapply(ls(),function(x){object.size(get(x))}))

import argparse
import numpy as np

############################## ARGUMENTS AND OPTIONS ################################
parser = argparse.ArgumentParser()

parser.add_argument("--verbosity", help="increase output verbosity")

parser.add_argument("--verbose", help="verbose or not (flag)", action="store_true")

parser.add_argument("--legend_file", type=str, help="the filled legend file", default='tmp/Wollstein.hg19.legend.csv')

parser.add_argument("--asd_file", type=str, help="the asd file", default='tmp/Wollstein.hg19.asd')

parser.add_argument("--mds_file", type=str, help="the plot in pdf", default='tmp/Wollstein.hg19.mds.pdf')

parser.add_argument("--xvector", help="dimension 1 (an integer)", type = int, default=1)
parser.add_argument("--yvector", help="dimension 2 (an integer)", type = int, default=2)


parser.add_argument("--option_mds", help="perform mds", action="store_true")
parser.add_argument("--option_pca", help="perform pca", action="store_true")

parser.add_argument("--option_summary", help="display summary", action="store_true")
parser.add_argument("--option_compact", help="compact graph", action="store_true")

parser.add_argument("--title_plot", type=str, help="the title of the plot", default='worldwide')

args = parser.parse_args()
verbose = args.verbose

legend_file = args.legend_file
asd_file = args.asd_file
mds_file = args.mds_file

xvector = args.xvector
yvector = args.yvector

option_mds = args.option_mds
option_pca = args.option_pca

option_summary = args.option_summary
option_compact = args.option_compact

title_plot = args.title_plot

dico_options = {"verbose": verbose,"option_compact":option_compact,"option_summary":option_summary,"option_mds":option_mds,"option_pca":option_pca,"xvector":xvector,"yvector":yvector,"mds_file":mds_file,"asd_file":asd_file,"legend_file":legend_file}

options = '''parse
legend_file <- %(legend_file)s;
asd_file <- %(asd_file)s;
mds_file <- %(mds_file)s;
xvector <- %(xvector)i;
yvector <- %(yvector)i;
option_mds <- %(option_mds)s;
option_pca <- %(option_pca)s;
option_summary <- %(option_summary)s;
option_compact <- %(option_compact)s;
PlotTitle <- "Worldwide" #used only if option_compact
verbose<-%(verbose)s;
'''
if verbose: 
    print(options%dico_options)
########################## read in data  ##############################

## Identify individuals removed in legend file.
def read_data_files(asd_file,legend_file):
    # check different table readers
    # numpy.loadtxt; a,b,c = np.loadtxt('data_table.txt', skiprows=1, unpack=True)
    # numpy.genfromtxt; np.genfromtxt('data_table3.txt', skip_header=1, missing_values=('MISSING','MISSING','MISSING'), filling_values=(-999,-999,-999))
    # asciitable; >>> np.genfromtxt('data_table3.txt', skip_header=1, missing_values=('MISSING','MISSING','MISSING'), filling_values=(-999,-999,-999))
    # atpy; t = atpy.Table('data_table.txt', type='ascii'); t_new = t.where((t['Pl. Mass'] > 5) & (t['Pl. Period'] < 100))
    asd_orig = np.loadtxt(asd_file);
    print(asd_orig)
    print
    return None

## Identify individuals with no shared markers
def identify_individuals_removed_in_legend_file():
    return None

## Unpack data from asd_orig
def unpack_data_from_asd_orig():
    return None

## Remove individuals
def remove_individuals():
    return None

## Population legends is all legends except individual and empty lines
def compute_legend_pop():
    return None
def compute_sample_reference_list():
    return None

## convert the "B" into ascii(B) ----> check, must be the R 
## Perform scaling
def perform_scaling():
    return None

## Determine the output file format and open accordingly
def output_format(): #open_plot_file
    return None

## Plot legend
def plot_legend():
    return None

## Plot populations
def plot_populations(): # decide in there if MDS or PCA
    return None

## Plot individuals
def plot_individuals(): # only if not compact requested
    return None

## Compute summary information
def summary_information():
    return None

## run all?
def all():
    read_data_files(asd_file,legend_file)
    return None

all()
