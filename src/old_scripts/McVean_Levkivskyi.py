#!/usr/bin/env python

from __future__ import division

import simul_ms
import numpy as np
import pylab as py 

def lambda_McVean(T,D,n1,n2):
    n = n1+n2
    return (1+2*D*n1*n2/n)/T

def lambda_Levkivskyi(T,D,n1,n2):
    n=n1+n2
    return (2+4*n1*n2/n*D*(D+2))/T**2

D = 0.1
n1 = 40
n2 = 60



lM = []
lL = []
ratio = []
Ds = np.linspace(0,3,20)
for dd in Ds:

    params,data,tree_lengths = simul_ms.ms_two_pops(n1=n1, n2=n2, D=dd, nreps=1000)

    T = np.average(tree_lengths)

    Mc = lambda_McVean(T,D,n1,n2)
    Lev = lambda_Levkivskyi(T,D,n1,n2)

    lM.append(Mc)
    lL.append(Lev)
    
    ratio.append(Mc/Lev)

py.title("2 populations, first eigenvalue")

py.plot(Ds,lM,'o',color='orange',label='McVean')

py.plot(Ds,lL,'v',color='blue',label='ours')

py.plot(Ds,ratio,'x',color='black',label='McVean/ours')

py.xlabel("D")

py.legend(loc='best')
py.savefig("lambda_McVean_ours_0_3.pdf")
py.show()



