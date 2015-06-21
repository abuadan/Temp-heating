from __future__ import division 
import matplotlib.pylab as plt
import sympy as sy
from sympy import besselk, log, nsolve
import math
import numpy as np
import scipy.special as sp
#These are Global vaiables.
Tc = 9.25
Tb = 4.2 #Tb represents the bath temperature 
tb = Tb / Tc
print 'Value of tb is:'+' '+ str(tb)
L = 100*10**-9
W = 80*10**-9
n = 2*10**-6
r1 = W / 2
print 'value of r1 is:'+' '+ str(r1)
x1 = r1 /n
print 'value of x1 is'+ ' '+ str(x1)
#S and B are varibales from the Hazara paper on the right hand coloume of page 2 
S = math.sqrt(3 / (1+ tb + tb**2))
print 'value of S is:'+' ' + str(S)
B = (math.pi / 2)*(1+(L/W))
print 'value of B is:'+' ' + str(B)
# The follwoing is done to calculate x0 from equation 13 in page 4 from the Hazara paper since T=Tc at x0 and therefore t = 1 rearraging equation 13
#gives the new equation t 
guess = [x /10 for x in range (1, 2, 1)] #The guess needs to be a small value to ensure correct answer, since we Know the value of x1 then it is safe to assume a an 
#initial guess near that value, the important thing to note is the answer for x0 must be larger than x1. when the guess range is changed it alteres the value of t and X0
x0 = sy.symbols('x0')
i = (S*x0)*(1-tb**3)*besselk(0, x0) /  3*(log(x0/x1)+B) * besselk(0, S*x0)
t = - i*( (log(x0/x1)*B)**2 - (log(x0/x1)+B)**2 )
X0 = nsolve(t, guess)
print 'value of x0 is:' +' '+ str(X0)
#Once the value of x0 is known the first calculations look into the normal region of the Weak link the region that is x1<= x < x0
#First a few variables are declared a few variables and list
x=x1
Temp=[]
X=[]
while x<X0:
    x +=0.0001
    print 'value of x is:' +' '+ str(x)
    E = ((S*float(X0))*(1-tb**3) *sp.k1(float(X0)*S)) / ( 3*(np.log(float(X0)/x1)+B) * sp.k0(S*float(X0)))
    print 'value of i is:' +' '+ str(E)
    t = math.sqrt(1 - E * ( ((np.log(x/x1)+B)**2) - ((np.log(float(X0)/x1)+B)**2) ))
    print t
    T = t*Tc
    Temp.append(T)
    X.append(x)
    print t
    print 'The Temperature at x='+' '+str(x)+' is '+str(T)
    if x==X0:
        break
# The same process is repeated for the superconducting region thats everywhere beyond x0 where the temperature goes to the bath temperature 
x2=X0
XX = 40*X0
Temp1=[]
x00=[]
while x2<XX:
    x2+=0.0001  
    t1= tb**3 + ((1-tb**3) / sp.k0(float(X0)*S))*sp.k0(float(x2)*S)
    T1= t1**(1/3)*Tc
    Temp.append(T1)
    X.append(x2)
    print 'The Temperature at x='+' '+str(x2)+' is '+str(T1)
    if x2==XX:
        break
plt.plot(X, Temp,'b-', linewidth=1.75)
plt.plot(X[0], Temp[0], 'ro')
plt.plot(X[444], Temp[444], 'go')
plt.ylabel('Temperature / K', fontsize=16)
plt.xlabel('x', fontsize=16)

plt.show()