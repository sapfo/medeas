#!/usr/bin/env python

"""
Created on Tue Oct  6 18:56:46 CEST 2015

@author: sapfo
"""
from __future__ import division

import numpy as np
import matplotlib.pyplot as pl
from textwrap import dedent
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.legend_handler import HandlerLine2D
import random

T = 2
RANDOM = 3

TEST = '4pops-grad'
THEORY = True

# --------------- setup ---------------

PRESETS = {

'2pops'     : ([20, 40],             # population sizes
               ['b', 'g'],           # colors for plotting
               '''
               d[0, 1] = 8
               '''),                 # distances

'3pops'     : ([50, 30, 20],
               ['b', 'g', 'r'],
               '''
               d[0, 1] = 3
               d[0, 2] = d[1, 2] = 8
               '''),

'4pops-pair' : ([110, 20, 30, 40],
                ['b', 'g', 'r', 'c'],
                '''
                d[0, 1] = 3
                d[2, 3] = 5
                d[0, 2] = d[0, 3] = d[1, 2] = d[1, 3] = 8
                '''),

'4pops-grad' : ([110, 20, 30, 40],
                ['b', 'g', 'r', 'c'],
                '''
                d[0, 1] = 3
                d[0, 2] = d[1, 2] = 5
                d[0, 3] = d[1, 3] = d[2, 3] = 8
                '''),

}

ns, cs, ds = PRESETS[TEST]
N = len(ns)                        # number of pops
n = sum(ns)                        # total number of individs

d = np.zeros((N, N))
exec(dedent(ds), globals())

for p in range(N):
    for q in range(p, N):
        d[q, p] = d[p, q]            # distances must be symmetric

# --------------- theoretical results ---------------

def two_pops(n1=5, n2=2, D=0.1, T=3):
    '''
    Calculate the largest eigenvector for two populations
    '''
    params = dict(zip(['n1','n2','D','T'],[n1,n2,D,T]))

    n=1.*(n1+n2)

    lambda_1 = (2/n/T**2)*(2*n1*n2*D*(D+2) + n)
    lambda_2 = (2/T**2)
    eigenvalues = [lambda_1]+int((n-2))*[lambda_2]
	
    x = np.sqrt(n2/n/n1)
    y = np.sqrt(n1/n/n2)
     
    v1 = n1*[np.sqrt(n2/n/n1)]+n2*[np.sqrt(n1/n/n2)]	
     
    eigenvectors = [v1]
    return params,eigenvalues, eigenvectors
 	

#    return [-np.sqrt(lambda_1)*x, np.sqrt(lambda_1)*y]


def three_pops(n1=2, n2=3, n3=3, D=0.1, D1=0.2, T=0.3):

    '''
    Calculate two largest eigenvectors for three populations

    ASSUME n1!=n2
    '''
    params = dict(zip(['n1','n2','n3','D','D1','T'],[n1,n2,n3,D,D1,T]))

    if n1==n2:
        raise ValueError('We have assumed n1 != n2, here n1=%i, n2=%i'%(n1,n2))
    n=1.*(n1+n2+n3); n1 = 1.*n1; n2 = 1.*n2; n3 = 1.*n3
    T = 1.*(T)
    D = 1.*D; D1= 1.*D1

    dt1 = D1*(D1+2)
    dt = D*(D+2)
    
    R_sqrt = np.sqrt(dt**2    * (n1+n2)**2*n3**2
               +dt1**2   * n1*n2*(n1+n3)*(n2+n3)
               -2*dt*dt1 * n1*n2*n3*(n+n3))
    	
    lambdat_p = -(dt1*n1*n2 + dt*(n1+n2)*n3 + R_sqrt)/n
    lambdat_m = -(dt1*n1*n2 + dt*(n1+n2)*n3 - R_sqrt)/n
    
    lambda_p = (2/T**2)*(1 - lambdat_p)
    lambda_m = (2/T**2)*(1 - lambdat_m)
	
    lambda_3 = (2/T**2)
    #print T
    eigenvalues = [lambda_p,lambda_m]+int((n-3))*[lambda_3]

    xp = -(-dt*(n1+n2)*n3 + dt1*(n2+n3)*n1 + R_sqrt)/(dt1*(n1-n2)*n1)
    yp =  (-dt*(n1+n2)*n3 + dt1*(n1+n3)*n2 + R_sqrt)/(dt1*(n1-n2)*n2)
    xm =  ( dt*(n1+n2)*n3 - dt1*(n2+n3)*n1 + R_sqrt)/(dt1*(n1-n2)*n1)
    ym = -( dt*(n1+n2)*n3 - dt1*(n1+n3)*n2 + R_sqrt)/(dt1*(n1-n2)*n2)
    
    print (dt1*(n1-n2)*n1)

    vp = int(n1)*[xp]+int(n2)*[yp]+int(n3)*[1]	
    vm = int(n1)*[xm]+int(n2)*[ym]+int(n3)*[1]
		
    eigenvectors = [vp,vm] 
    return params, eigenvalues, eigenvectors

def T_D_two_pops(eigenvalues=[1,0.1,0.1,0],n1=20,n2=20,gen_time=29,Ne=10000,diploid=2,verbose=True):

    if verbose: print "eigenvalues: ",eigenvalues

    lambda_1 = eigenvalues[0]
    lambda_av = np.average(eigenvalues[1:-1])

    if verbose:
        print "lambda_1: ",lambda_1
        print "lambda_av: ",lambda_av

    T = np.sqrt(2./lambda_av)
    
    a = 1; b = 2; c = -((lambda_1/lambda_av)-1)*((n1+n2)/(2*n1*n2)) 

    D = (-b+np.sqrt(b**2-(4*a*c)))/(2*a)

    Drescaled = 1.*D*gen_time*diploid*Ne

    if verbose: 
        print "T: ",T
        print "D: ",D
        print "Drescaled(years): ",Drescaled    
    return T,D,Drescaled

def quadratic(a,b,c):
    pos = (-b+np.sqrt(b**2-(4*a*c)))/(2*a)
    neg = (-b+np.sqrt(b**2-(4*a*c)))/(2*a)
    return pos,neg
    
def T_D_D1_three_pops(eigenvalues=[1,0.1,0.1,0],n1=20,n2=20,n3=20,gen_time=29,Ne=10000,diploid=2,verbose=True):
    if verbose: 
        print "eigenvalues: ",eigenvalues
        print "n1, n2, n3: ",n1,n2,n3
    lambda_p = eigenvalues[0]
    lambda_m = eigenvalues[1]
    lambda_av = np.average(eigenvalues[2:-1])

    if verbose:
        print "lambda_1: ",lambda_p
        print "lambda_2: ",lambda_m
        print "lambda_av: ",lambda_av
    
    T = np.sqrt(2./lambda_av)
               
    if n1!=n2 or n1!=n3 or n2!=n3:
        raise ValueError("this has only been worked out for n1==n2==n3\n YOU CANT RUN IT IF NOT")
    else:
        n = n1+n2+n3
        N = n1
        
        d = (-8+(lambda_m+3*lambda_p)*T**2)/(8*N)
        d1 = (-2+(lambda_m*T**2))/(2*N)
    
        if verbose:
            print "d: ",d
            print "d1: ",d1
            
        Dpos,Dmin = quadratic(1,2,-d)  
        D1pos,D1min = quadratic(1,2,-d1)
        
        D = max(Dpos,Dmin)         
        D1 = max(D1pos,D1min)            
                
        Drescaled = 1.*D*gen_time*diploid*Ne
        D1rescaled = 1.*D1*gen_time*diploid*Ne
        
        if verbose:     
            print "T: ",T
            
            print "D: ",D
            print "D1: ",D1
            
            print "Drescaled(years): ",Drescaled    
            print "D1rescaled(years): ",D1rescaled    
            
        return T,D,D1,Drescaled,D1rescaled

def dist_two_populations(T=0.1, D=0.3, n1=2, n2=2):
    n = n1+n2
    dist_squared = 2./(T**2) * (2*D*(D+2)+n/(n1*n2))
    return  np.sqrt(dist_squared)

def dist_two_populations_resc_diag(T=0.1, D=0.3, n1=2, n2=2):
    n = n1+n2
    dist_squared = 2./(T**2) * (2*D*(D+2)+(2./n)+1)
    return  np.sqrt(dist_squared)

def dist_two_populations_resc_all(T=0.1, D=0.3, n1=2, n2=2):
    n = n1+n2
    dist_squared = 2./(T**2) * (D*(D+2)*(1-2/n))
    return  np.sqrt(dist_squared)




#T_D_D1_three_pops(eigenvalues = [0.8,0.2,0.02,0], n1=20,n2=20,n3=20)


#params,eigenval,eigenvec = two_pops(n1=10, n2=10, D=0.1, T=3)
#params,eigenval,eigenvec = three_pops()

#print eigenval
#print len(eigenval)
#print eigenvec
#print params

