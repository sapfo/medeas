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
parser = argparse.ArgumentParser()

parser.add_argument("--verbosity", help="increase output verbosity")
parser.add_argument("--verbose", help="verbose or not (flag)", action="store_true")

parser.add_argument("--legend_file", type=str, help="the filled legend file", default='../doc/HDGP.Reference.legend.filled.csv')

parser.add_argument("--legend_file", type=str, help="the asd file", default='../doc/HDGP.Reference.asd')

args = parser.parse_args()

verbose = args.verbose

legend_file = args.legend_file
asd_file = args.legend_file

b = '''parse

legend_file <- "tmp/Han_700000reads_hg19,Wollstein.hg19,.legend.filled.csv";
asd_file <- "tmp/Han_700000reads_hg19,Wollstein.hg19,.asd";
mds_file <- "/tmp/mike.pdf"
xvector <- 1;
yvector <- 2;
option_mds <- T;
option_pca <- T;
option_summary <- T;
option_compact <- F;
PlotTitle <- "Worldwide" #used only if option_compact
'''
if verbose: 
    print(b)
