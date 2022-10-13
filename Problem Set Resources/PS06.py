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

#the integrand of the integrated Friedmann Equation
def Ht(a,omega_0,omega_r,omega_m,omega_lam):
    H = omega_r/(a*a)+omega_m/a+omega_lam*a*a+(1.0-omega_0)
    H = np.sqrt(H)
    H = 1/H
    return H

#define the speed of light (Mpc/year)
c = 3.066e-7

#define the benchmark universe
omega_0 = 1.0
omega_m = 0.306
omega_lam = 0.692
omega_r = 9.03e-5

#inverse Hubble parameter H0 in years
H0inv = 1.4432e10

#calculate the age of the Universe
#integrate from a = 0 to a=1
Ht0 = quad(Ht,0,1.0,args=(omega_0,omega_r,omega_m,omega_lam))[0]

t0 = Ht0*H0inv

#answer to 1c
print('age of the Universe', t0)
#answer to 2c
print( 'age at decoupling', H0inv*quad(Ht,0,1.0/1101,args=(omega_0,omega_r,omega_m,omega_lam))[0])

#create a logrithmically populated list for the upper limits of scale factor integration
a_space = np.logspace(-6,9,10000)

#create an empty list that will be populated by the times for each integration
t_space = np.array([])

#perform a loop where I integrate from a=0 to a for each a in the list
#then add the result to tlist
for a in a_space:
    Ht_a = quad(Ht,0,a,args=(omega_0,omega_r,omega_m,omega_lam))[0]
    t_space = np.concatenate([t_space,[H0inv*Ht_a]])

#make a plot of tlist versus alist
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

plt.xlabel('$t$ (years)',fontsize=20)
plt.ylabel(r'$a$',fontsize=20)

#set the axis limits in years
ax.set_xlim([1e3,4e11])

#log log plot
ax.set_xscale('log')
ax.set_yscale('log')

ax.axhline(y=1,color='k',linestyle=':')
ax.axvline(x=t0,color='k',linestyle=':')
plt.plot(t_space,a_space,'b',linewidth=1.5,label=r'benchmark Universe')

plt.legend(loc='upper left')
plt.savefig("PS06_b.pdf",format='pdf')


#make a plot of tlist versus alist
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

plt.xlabel('$t$ (years)',fontsize=20)
plt.ylabel(r'$z$',fontsize=20)

#set the axis limits in years
ax.set_xlim([1e3,4e11])

#log log plot
ax.set_xscale('log')
ax.set_yscale('log')

#ax.axhline(y=1,color='k',linestyle=':')
ax.axvline(x=t0,color='k',linestyle=':')
plt.plot(t_space,1/a_space-1,'b',linewidth=1.5,label=r'benchmark Universe')

plt.legend(loc='upper left')
plt.savefig("PS06_b_z.pdf",format='pdf')

#now calculate the age of the Universe without radiation
Ht0_norad = quad(Ht,0,1.0,args=(omega_0,0,omega_m,omega_lam))[0]

t0_norad = Ht0_norad*H0inv

print('age of the Universe (no radiation)', t0_norad)
#take the different in ages
#answer to 1d
print('t0-t0_norad', t0-t0_norad)

#define an interpolating function of a as a function of time
def a_interp(t):

    return np.interp(t,t_space,a_space)

#conversion between inverse years and km/s/Mpc
yr_to_kmsMpc = 9.78e11


H_space = np.array([])
for t in t_space:
#calculate the derivative of a_interp, divide by a
    h = yr_to_kmsMpc*derivative(a_interp,t)/a_interp(t)
    H_space = np.concatenate([H_space,[h]])

#make a plot of tlist versus Hlist
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

#set the axis limits in years
ax.set_xlim([1e3,4e11])

plt.xlabel('$t$ (years)',fontsize=20)
plt.ylabel(r'$H$ (km/s/Mpc)',fontsize=20)

#log log plot
ax.set_xscale('log')
ax.set_yscale('log')

ax.axhline(y=1,color='k',linestyle=':')
ax.axvline(x=t0,color='k',linestyle=':')
plt.plot(t_space,H_space,'b',linewidth=1.5,label=r'benchmark Universe')

plt.legend(loc='upper left')
plt.savefig("PS06_e.pdf",format='pdf')


#define an interpolating function of 1/a as a function of time
def ainv_interp(t):

    return 1.0/np.interp(t,t_space,a_space)

dp_space = np.array([])

#integrate c*1/a from t to t0
for t in t_space:
    d = quad(ainv_interp,t,t0)[0]
    dp_space = np.concatenate([dp_space,[c*d]])


fig = plt.figure()
ax = fig.add_subplot(1,1,1)

plt.xlabel('$z$',fontsize=20)
plt.ylabel(r'$d_p(z)$ (Mpc)',fontsize=20)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-3,1e4])

plt.plot(1/a_space-1,dp_space,'b',label='$d_p(t_0,t_e)$')

plt.legend(loc='lower left')
plt.savefig("PS06_f.pdf",format='pdf')


fig = plt.figure()
ax = fig.add_subplot(1,1,1)

plt.xlabel('$z$',fontsize=20)
plt.ylabel(r'$d_p(z)$ (Mpc)',fontsize=20)

ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xlim([1e-3,1e4])

plt.plot(1/a_space-1,dp_space*a_space,'b',label='$d_p(t_e)$')

plt.legend(loc='lower left')
plt.savefig("PS06_g.pdf",format='pdf')


#calculate the age of the Universe
#when z = 1100
Ht_cmb = quad(Ht,0,1.0/1101,args=(omega_0,omega_r,omega_m,omega_lam))[0]

t_cmb = Ht_cmb*H0inv


Ht_rm = quad(Ht,0,1.0/3500,args=(omega_0,omega_r,omega_m,omega_lam))[0]

t_rm = Ht_rm*H0inv

print(t_rm)

#answers to 2c/d/e
print('age of the Universe at CMB decoupling:', t_cmb)
print('comoving distance to CMB surface:', c*quad(ainv_interp,t_cmb,t0)[0])
print('distance to CMB surface at emission time:', c*quad(ainv_interp,t_cmb,t0)[0]/1101)

print('age of 35', H0inv*quad(Ht,0,6.8e-6,args=(omega_0,omega_r,omega_m,omega_lam))[0])
print('distance to CMB surface at emission time:', c*quad(ainv_interp,35.0,t0)[0]*6.8e-6)
plt.show()
