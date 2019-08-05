# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 10:47:39 2016

@author: ben
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
"""
# Thermal Constants (Van der Veen 2013 - from Yen CRREL 1981)
spy  = 31556926.            #Seconds per year
rhoi = 917.                 #Bulk density of ice (kg/m3)
Kice = 1.09e-6              #Cold ice diffusivity (m2/sec)
KiYr = Kice*spy             #Ice Diffusivity (m2/yr)

### Draw a springtime temperature profile based on the equation in Aschwandens 2012 Thermodynmaics course
w = 2*np.pi                 # period of surface oscillation
ds = np.arange(0,20.1,0.1)  # depth below the surface (m)
MoY = 2.                    # Month of year for the desired temperature profile
Tc = 263.39                  # Temperature profile oscillates around the convergence temp T0 (K)
Ta = 10.0                   # Size of annual oscillations at the surface (K)

def Temp(d):
    return Tc + Ta*np.exp(-d*np.sqrt(w/(2.*KiYr)))*np.sin((w*(MoY/12.)) - d*np.sqrt(w/(2.*KiYr)))

### Calculate the rate factor , A, over a range of temperatures based on the equation in Van der Veen 2013 (p 33)
A0 = 9.302e7               # rate factor constant (kPa^-3 yr^-1)
Q = 78800.                    # Activation energy for creep (J/mol)
R = 8.321                   # Gas Constant (J/mol K)                     
Tr = 273.39                 # Reference Temperature (K)
K = 1.17                    # Power constant (unitless)
C = 0.16612                 # constant (K**k)

def RateFactor(T):
    return A0*np.exp(-(Q/(R*T))+(3*C/((Tr-T)**K)))  # rate factor (kPa^-3 yr^-1)  

######################################################################################
### Fracture Mechanics  

# Calculate the effective stress based on some constant effective strain rate
g = 9.81                    #Gravity (m/s2)
K1c = 100.                   # Fracture Toughness (kPa m^.5) Sculson and Duval 2009
edot = 10e-9                 # effective strain rate (s-1, see Sculson and Duval pg 79 low strain rate 10e-7, high 10e-3)
n = 3.                      # Power constant in Glen's Law

# Determine the total applied stress at some depth, d, based on the constant strain rate
def AppStress(d):
    Rxx = ((edot*spy)/RateFactor(Temp(d)))**(1/n)
    Lxx = -rhoi*g*d/1000.
    return Rxx + Lxx

# Geometry-dependent constant, should be of order unity
# Penny Crack (Lawn 1993)
psi = 2./(np.pi)**.5
# line crack below surface (Tada 1973)
def Gfactor(c,d):
    return 2/np.sqrt(np.pi)#np.sqrt(1/(np.cos((np.pi*c)/(2*d))))
# Fracture Toughness of Ice
K1c = 100.

# Determine the stress intensity factor by assuming that the stress over the crack is all the same
def K1point(c,d):
    return Gfactor(c,d)*AppStress(d)*c**0.5

# Determine the stress intensity factor by integrating the applied stress over the depth of the crack (Tada 1973 p. 2.27)
def f(d):
    return Gfactor(c,d)*AppStress(d)
def K1top(c,d):
    return 2/(np.sqrt(np.pi*c))*(quad(f,d-c,d)[0])
def K1bottom(c,d):
    return 2/(np.sqrt(np.pi*c))*(quad(f,d,d+c)[0])
    

# crack sizes = 2c (m)
cs = np.arange(0,50.1,0.1) 
cs = np.insert(cs,1,0.01)           
# calculate stress intensity factors
SIFtop = np.array([[K1top(c,d) for d in ds] for c in cs])
SIFbottom = np.array([[K1bottom(c,d) for d in ds] for c in cs])
SIFdiff = SIFbottom-SIFtop
# create a grid for contouring
DS,CS = np.meshgrid(ds,cs)
"""
fig = plt.figure(figsize=(9,6))

ax1 = plt.subplot(121)
plt.title('top of crack')
plt.xlabel('Crack Radius ($m$)')
plt.ylabel('Depth Below Surface ($m$)')
plt.ylim(max(ds),min(ds))
plt.xlim(0.01,1)#max(cs))
plt.contourf(CS,DS,SIFtop,np.linspace(np.nanmin(SIFbottom),np.nanmax(SIFbottom),100),cmap=plt.get_cmap('RdYlBu_r'))
#plt.colorbar(ticks=np.arange(0,np.nanmax(SIFtop),25),label="Stress Intensity Factor $kPa * m^{1/2}$")
plt.hold(True)
plt.contour(CS,DS,SIFtop,[K1c],colors='k',linewidths=2)

ax2 = plt.subplot(122,sharey=ax1)
plt.title('bottom of crack')
plt.xlabel('Crack Radius ($m$)')
#plt.ylabel('Depth Below Surface ($m$)')
plt.ylim(max(ds),min(ds))
plt.xlim(0.01,1)#max(cs))
plt.contourf(CS,DS,SIFbottom,np.linspace(np.nanmin(SIFbottom),np.nanmax(SIFbottom),100),cmap=plt.get_cmap('RdYlBu_r'))
plt.colorbar(ticks=np.arange(np.nanmin(SIFbottom),np.nanmax(SIFbottom),100),label="Stress Intensity Factor $kPa * m^{1/2}$")
plt.hold(True)
plt.contour(CS,DS,SIFbottom,[K1c],colors='k',linewidths=2)

plt.savefig('TopBottom_2.png')
#"""