#!/usr/bin/env python


from __future__ import division
 
import numpy as np
import pylab

def cmdscale(D):
    """                                                                                       
    Classical multidimensional scaling (MDS)                                                  
                                                                                               
    Parameters                                                                                
    ----------                                                                                
    D : (n, n) array                                                                          
        Symmetric distance matrix.                                                            
                                                                                               
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
    # Number of points                                                                        
    n = len(D)
 
    # Centering matrix                                                                        
    H = np.eye(n) - np.ones((n, n))/n
 
    # YY^T                                                                                    
    B = -H.dot(D**2).dot(H)/2
 
    # Diagonalize                                                                             
    evals, evecs = np.linalg.eigh(B)
 
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

if 0:
    T = 2
    Delta = 2/(T)*np.array([[0,1,1],[1,0,1],[1,1,0]])
    Delta = np.array([[0,1,1],[1,0,1],[1,1,0]])


    evals, evecs, Y = cmdscale(Delta)
    print "evals ", evals
    print "coords ",Y
    
    print "coords ",Y[0,:]
    pylab.plot(Y[:,0],Y[:,1],'o')
    pylab.show()

