# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:13:10 2016

@author: ben

This script recreates what Nath and Vaughan (2003) worked up as theory for 
subsurface crevasses in firn. 
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Constants (Van der Veen 2013 - from Yen CRREL 1981)
g = 9.81                    #Gravity (m/s2)
rhoi = 917.                 #Bulk density of ice (kg/m3)
K1c = 100                   # Fracture Toughness (kPa m^.5) Sculson and Duval 2009
rhos = 450.                 # Surface Snow Density
ds = np.arange(50.)

# define a depth dependent density of firn
def rho(depth):
    return rhoi - (rhoi-rhos)*np.exp(-C*depth)
# calculate the lithostatic stress (kPa) based on the density profile where C is the densification rate
L = np.array([[-(g/1000.)*quad(rho,0,d)[0] for d in ds] for C in [.015,.03]])

# Define fracture toughness as a function of density (from Rist et al. 1999)
# x1 and x2 are constants, x2 changes for upper and lower bounds
x1 = .2567
K1c = np.array([[x1*rho(d)+x2 for d in ds] for x2 in [-80.7-18.7,-80.7+18.7]])




### Nath and Vaughan Scenario #1, Constant Dynamic Tensile Stress with Depth
"""
beta = 1.12                 # Geometric constant (for now this is a guess from Van der Veen because I couldn't get their reference (Isida 1966))
sigma = np.array([Rxx+L for Rxx in [100,250,400]])              # Applied stress (kPa)
sigma[sigma<0.0]=0.0


a = [(K1c[0]/(beta*sig[1]*np.sqrt(np.pi)))**2 for sig in sigma]

plt.xlabel('Starter Crack Length ($m$)')
plt.ylabel('Depth Below Surface ($m$)')
plt.title('Dynamic Tensile Stresses')
plt.plot(a[0],ds,'k',label='100 kPa')
plt.plot(a[1],ds,'k--',label='250 kPa')
plt.plot(a[2],ds,'k:',label='400 kPa')

plt.xlim(0,15)
plt.ylim(50,0)
plt.grid('on')
plt.legend()
plt.savefig('NandV_7c.png')#"""


### Nath and Vaughan Scenario #2, Constant Strain Rate with Depth

ABs = np.array([[[.455,18.981e-3],[6.4e3,3.054e-3]],[[1.929,18.002e-3],[1.285e4,3.243e-3]],
       [[.15,23.021e-3],[2.316e4,3.144e-3]]])

def Rxx(rho,A,B):
    return A*np.exp(B*rho)

C = 0.03
L = np.array([-(g/1000.)*quad(rho,0,d)[0] for d in ds])
R = np.zeros((3,len(ds)))
for i in range(3):
    for d in ds:
        if rho(d)<600.:
            R[i,d] = Rxx(rho(d),ABs[i,0,0],ABs[i,0,1])/1000.
        else:
            R[i,d] = Rxx(rho(d),ABs[i,1,0],ABs[i,1,1])/1000.

"""
fig = plt.figure()
plt.title('Density Dependent Stresses')
plt.xlim(-300,300)
plt.xlabel('Stress ($kPa$)')
plt.ylim(60,0)
plt.ylabel('Depth Below Surface ($m$)')
plt.plot(L,ds,'b-',label='Lithostatic')
plt.plot(R[0],ds,'k-',label='100 kPa')
plt.plot(R[1],ds,'k-',label='250 kPa')
plt.plot(R[2],ds,'k-',label='400 kPa')
plt.grid('on')
plt.legend(loc=3)
plt.savefig('NandV_8.png')"""

beta = 1.12                 # Geometric constant (for now this is a guess from Van der Veen because I couldn't get their reference (Isida 1966))
sigma = R+L
sigma[sigma<0.0]=0.0


a = [(K1c[0]/(beta*sig*np.sqrt(np.pi)))**2 for sig in sigma]

plt.xlabel('Starter Crack Length ($m$)')
plt.ylabel('Depth Below Surface ($m$)')
plt.title('Dynamic Tensile Stresses')
plt.plot(a[0],ds,'k--',label='100 kPa')
plt.plot(a[1],ds,'k--',label='250 kPa')
plt.plot(a[2],ds,'k:',label='400 kPa')

plt.xlim(0,15)
plt.ylim(50,0)
plt.grid('on')
plt.legend(loc=4)
plt.savefig('NandV_9c.png')
