# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 23:14:36 2015

@author: ivan
"""
from __future__ import division

from sympy import *

n1,n2,n3,d,d1 = symbols('n1 n2 n3 d d1')

n = n1+n2+n3
att = (n1*(n1-1)+n2*(n2-1)+n3*(n3-1)+2*n1*n2*d1+2*(n1+n2)*n3*d)/n**2

a1 = (n1-1+n2*d1+n3*d)/n
a2 = (n1*d1-1+n2+n3*d)/n
a3 = (n1*d-1+n2*d+n3)/n

a = att-2*a1
b = att-2*a2
c = att-2*a3
f = d1-a1-a2+att
e = d-a1-a3+att
g = d-a2-a3+att

aa = n1*(a+1)
bb = n2*(b+1)
cc = n3*(c+1)

M = Matrix([[aa,n2*f,n3*e],[n1*f,bb,n3*g],[n1*e,n2*g,cc]])
M = simplify(M)
print(M)
v = M.eigenvects()