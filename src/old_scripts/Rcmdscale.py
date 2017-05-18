#!/usr/bin/env python

"""
Created Wed Oct  7 14:56:43 CEST 2015

@author: sapfo
"""

"""
points: a matrix with up to k columns whose rows give the
          coordinates of the points chosen to represent the
          dissimilarities.

eig: the n eigenvalues computed during the scaling process if
       eig is true.  *NB*: versions of R before 2.12.1 returned
       only k but were documented to return n - 1.

x: the doubly centered distance matrix if x.ret is true....that must be the eigenvectors

"""
import numpy as np
import rpy2.robjects as ro
import sys

r = ro.r

# read from csv file

Delta = np.array([[0,1,2],[1,0,3],[2,3,0]])

dataf = r('read.table("table_essai.txt")')

B = r.cmdscale(dataf)

print B

#cmdscale(read.table("table_essai.txt"),k=2,eig=T,add=F,x.ret=T)$eig
#cmdscale(read.table("table_essai.txt"),k=2,eig=T,add=F,x.ret=T)$points
#cmdscale(read.table("table_essai.txt"),k=2,eig=T,add=F,x.ret=T)$x


#D: divergence time
#Delta: dissimilarity matrix
#total population size







sys.exit()

'''
n=5

D=3

T=4.5


Delta <- 2/T*matrix(c(0,1,1,(1+D),(1+D),1,0,1,1+D,1+D,1,1,0,1+D,1+D,1+D,1+D,1+D,0,1,1+D,1+D,1+D,1,0),n,n,byrow=T)


A=-1/2*(Delta^2)

#B=A-(2*1/n*(sum(A[1,])))+(1/n)^2*sum(A) WRONG...if n_a!=n_b

lambdas = eigen(B)$values

cmd<-cmdscale(Delta,k=n-1,eig=T,add=F,x.ret=T)

loc = cmd$points

x <- loc[, 1]
y <- loc[, 2]
z <- loc[, 3]
plot(x, y,  xlab = "x", ylab = "y", asp = 1, axes = TRUE, main = "pop subdivision dim1, dim2")

plot(x, z,  xlab = "x", ylab = "z", asp = 1, axes = TRUE, main = "pop subdivision")
plot(y, z,  xlab = "y", ylab = "z", asp = 1, axes = TRUE, main = "pop subdivision")



eigenvalue = -0.547723
eigenvector = c(0.365148,...,-0.547723)
distance between the extreme points is: 1.74507
'''
