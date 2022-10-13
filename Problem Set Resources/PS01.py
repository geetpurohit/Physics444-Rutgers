#!/usr/bin/env python

import math
import fileinput
import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.optimize
import math
import matplotlib.pyplot as plt
from scipy.integrate import quad, dblquad, tplquad, nquad
from scipy.misc import derivative

from matplotlib import rcParams

#this makes plots look prettier
rcParams.update({'figure.autolayout': True})
plt.rc('font', family='serif')


#for part a I first define a function
def f(x):
    A = 0.8
    B = 0.2
    
    f = A*pow(x,2.0)+B*pow(x,4.0)
    f = pow(f,0.5)
    return f

#In python, I will need to define some array of a values over which I will eventually integrate
#the numpy linspace command does this easily
a_list = np.linspace(1,100,1001)

#now I'll create an empty list to fill with the integral
F_list = np.array([])

#integrating a function in python uses the "quad" command
for a in a_list:
    F = quad(f,0,a,args=())[0]
    #store the integral in the new list
    F_list = np.concatenate([F_list,[F]])

#now create a plot
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

plt.xlabel('$a$',fontsize=20)
plt.ylabel(r'$\int_0^a f(x) dx$',fontsize=20)

plt.plot(a_list,F_list,linewidth=1.5)

plt.savefig('PS01_1a.pdf',format='pdf')

#for part B, it is trivial to swap the x and y axes 
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

plt.ylabel('$a$',fontsize=20)
plt.xlabel(r'$\int_0^a f(x) dx$',fontsize=20)

plt.plot(F_list,a_list,linewidth=1.5)

plt.savefig('PS01_1b.pdf',format='pdf')

#For part C, there are several ways to proceed. I'm going to do something that is overkill, so you see how its done. 
#first I'll interpolate the function
def F_interp(a):

    return np.interp(a,a_list,F_list)

#then I'll calculate the derivative of this interpolated function
diff_list = np.array([])

for a in a_list:
    deriv = derivative(F_interp,a)
    diff_list = np.concatenate([diff_list,[deriv]])

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

plt.xlabel('$a$',fontsize=20)
plt.ylabel(r'$\frac{dF(a)}{da}$',fontsize=20)

plt.plot(a_list,diff_list,linewidth=1.5,label='numeric derivative')
plt.plot(a_list,f(a_list)+100,linewidth=1.5,color='r',label=r'f(a)+100')

plt.legend(loc='upper left')

plt.savefig('PS01_1c.pdf',format='pdf')

plt.show()
