from __future__ import division
import matplotlib.pylab as plt
import numpy as np
import scipy.special as sp
import scipy as sci
from scipy import integrate
import math
from sympy.functions import coth
Tc = 9.25
Tb = 7.2
t = Tb / Tc
Temperature = []
Temp1=[]
Temp0=[]
D=[]
d = []
D1=[]
d1 = []
n = 2*10**-6
L = 100*10**-9
W = 80*10**-9
a = 3*10**-2
s1 = W/ (2*n)
y1 = (L+(W/2)) / (2*n)
x0 = 0.015
r0 = 2*x0
s2 = r0 / n
y0 = (x0 / n)/1000000
print x0, y0, y1
A = ((W/n)**2) *(sp.kv(0, s1)+(math.pi / 2)*sp.kv(1,s1)*coth(y1))
B = ((W/n)**2) *(sp.iv(0, s1)+(math.pi / 2)*sp.iv(1,s1)*coth(y1))
print A, B
def t1(t):
    return  (t**-1)*sp.kv(0, s2)
def t2(t):
    return (t**-1)*sp.iv(0, s2)
print t2
temp = integrate.quad(t1, s1, s2)
Fk2 = (math.pi**-2) * np.array(temp[:0])
temp1 = integrate.quad(t2, s1, s2)
FI2 = (math.pi**-2) * np.array(temp1[:0])
print Fk2 , FI2
r1 = 0.0
while r1 < y1:
    C0 = sp.kv(0,s2)*(1 + (A*FI2)-(B*Fk2))/A
    #print C0
    D_ = 1 - B*Fk2 - A*Fk2*sp.iv(1, s1) / sp.kv(1, 1)
    #print D_
    r1 += 0.0001
    j = -1*r1
    D.append(r1)
    d.append(j)
    T = Tb + (Tc - Tb) * (sp.kv(0,s1) + (math.pi /2)* sp.kv(1, s1)*coth(r1))*(1- D_ * math.cosh(y1)) * (C0*A)
    Temp0.append(T)
    #print Temp0, r1