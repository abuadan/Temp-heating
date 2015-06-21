from __future__ import division
import matplotlib.pylab as plt
import numpy as np
import scipy.special as sp
import scipy as sci
import math
from sympy.functions import coth
Tc = 9.25
Tb = 4.2
Temperature = []
Temp1=[]
Temp0=[]
D=[]
d = []
D1=[]
d1 = []
n = 2*10**-6
L = 50*10**-6
W = 80*10**-9
a = 3*10**-2
s1 = W/ (2*n)
y1 = (L+(W/2)) / (2*n)
x0 = 0.015
y0 = (x0 / n)/1000000
print x0, y0, y1
g = (2/math.pi)*(sp.kv(0,s1))*(sp.kv(1,s1))**-1
G = (1+g*math.tanh(y1-y0))*(1+g*coth(y1-y0))**-1
print g, G
r1 = 0.0
while r1 < y0:
    r1 += 0.0001
    j = -1*r1
    D.append(r1)
    d.append(j)
    T = Tc + (Tc - Tb)*G*coth(y0)*coth(y1-y0)*(1-(math.cosh(r1)/math.cosh(y0)))
    Temp0.append(T)
    print T, r1
r2 = y0
while r2 < y1:
    r2 += 0.0001
    D1.append(r2)
    d1.append(-1*r2)
    H0 = math.sinh(y1-y0) + g*math.cosh(y1-y0)
    Hy = math.sinh(y1-r2) + g*math.cosh(y1-r2)
    T1 = Tb +(Tc-Tb)*Hy/H0   
    Temp1.append(T1)
    print T1, r2
Temperature = Temp0 + Temp1
positive = D + D1
negative = d + d1
plt.plot(positive, Temperature, 'b-')
plt.plot(D[74], Temp0[74], 'ro')
plt.show()
print len(Temp0)
print len(positive)
print len(negative)