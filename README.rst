======
medeas
======
listing of directories, 02.01.18
================================
main.py: main entry point
options.py
processing
__pycache__
README.rst
simulate.py: script that runs many other scripts
src
test

tracy.py: not using yet, it gives pvalues as a  function of sigma (for now using a fixed pvalue so there is a single value of the stats)

transcode.py: converts a filtered output of scrm (filtered for numbers only) to the format that medeas use and saves the seed to a hardcoded filename: _tmp_seed,txt
arguments: folder, whereiscrm, splittime /currently works with two pops, one D)

within transcode.py:
create fake_labs because medeas needs it

for now, n1, n2, theta: hard coded.
