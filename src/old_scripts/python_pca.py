#!/usr/bin/env python


from __future__ import division
 
import numpy as np
import pylab

def PCA(Z, verbose = 0):
    """                                                                                       
    PCA Mc Vean style                                                  
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    Z : (L, n) array                                                                          
        SNP array matrix. Matrix of 0s and 1s that represent the two alleles.                                                            
                                                                                               
    Returns                                                                                   
    -------                                                                                   
    Y : (n, p) array                                                                          
        Configuration matrix. Each column represents a dimension. Only the                    
        p dimensions corresponding to positive eigenvalues of B are returned.                 
        Note that each dimension is only determined up to an overall sign,                    
        corresponding to a reflection.                                                        
        dimension 0: Y[:,0]
        dimension 1: Y[:,1]
        dimension 2: Y[:,2]
        ...

                                                                                       
    e : (n,) array                                                                            
        Eigenvalues of B.                                                                     
                                                                                               
    """
    # Number of SNPs                                                                        
    L = len(Z)
    if verbose: print "number of SNPs: ",L
 
    # Centering matrix  remove the average of each row
    X = Z-np.average(Z,axis=1).reshape(L,1)
      

    # matrix to diagonalise     
    M = 1./L*(X.T.dot(X))      
                                                            
    # Diagonalize                                                                             
    evals, evecs = np.linalg.eigh(M)
 
    # Sort by eigenvalue in descending order                                                  
    idx   = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]
 
    # Compute the coordinates using positive-eigenvalued components only                      
    w, = np.where(evals > 0)
    L  = np.diag(np.sqrt(evals[w]))
    V  = evecs[:,w]
    Y  = V.dot(L)
 
    return evals,evecs,Y

'''
Data = np.array([[0,1],[1,0],[0,0]])
evals, evecs, Y = PCA(Data)
print "evals ", evals
print "coords ",Y

print "coords ",Y[0,:]
pylab.plot(Y[:,0],'o')
pylab.show()
'''
