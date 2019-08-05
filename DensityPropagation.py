# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 14:05:34 2016

@author: ben
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Constants (Van der Veen 2013 - from Yen CRREL 1981)
g = 9.81                    #Gravity (m/s2)
rhoi = 917.                 #Bulk density of ice (kg/m3)
K1c = 100                   # Fracture Toughness (kPa m^.5) Sculson and Duval 2009
rhos = 450.                 # Surface Snow Density
ds = np.arange(0,50.1,2)
C = 0.03

### Nath and Vaughan Scenario #2, Constant Strain Rate with Depth
beta = 1.12                 # Geometric constant (for now this is a guess from Van der Veen because I couldn't get their reference (Isida 1966))

# define a depth dependent density of firn
def rho(d):
    return rhoi - (rhoi-rhos)*np.exp(-C*d)

# define a depth dependent applied stress (tensile stress plus lithostatic stress)
def Tensile(d):
    if rho(d) < 600.:
        A, B = 0.15, 23.021e-3
    else:
        A, B = 2.316e4, 3.144e-3
    return A*np.exp(B*rho(d))/1000.
        
def Lithostatic(d):
    return -rho(d)*g*d/1000.

def sigma(d):
    stress = Tensile(d) + Lithostatic(d)
    if stress > 0.:
        return stress
    else:
        return 0.
        
# Determine the stress intensity factor by integrating the applied stress over the depth of the crack (Tada 1973 p. 2.27)
def Gfactor(c,d):
    return 1.12

def f(d):
    return Gfactor(c,d)*sigma(d)

def K1calc(c,d):
    return 2/(np.sqrt(np.pi*c))*(quad(f,d-c,d+c)[0])

def K1top(c,d):
    if c>d:
        return 0.
    else:
        return 2/(np.sqrt(np.pi*c))*(quad(f,d-c,d)[0])
    
def K1bottom(c,d):
    return 2/(np.sqrt(np.pi*c))*(quad(f,d,d+c)[0])
       
# Define fracture toughness as a function of density (from Rist et al. 1999)
# x1 and x2 are constants, x2 changes for upper and lower bounds
x1 = .2567
K1c = np.array([[x1*rho(d)+x2 for d in ds] for x2 in [-80.7-18.7,-80.7+18.7]])
K1c = 100.

# crack sizes = 2c (m)
cs = np.arange(0,50.1,2)  
cs = np.insert(cs,1,0.01)     
# calculate stress intensity factors
SIF = [[K1calc(c,d) for d in ds] for c in cs]
SIFtop = [[K1top(c,d) for d in ds] for c in cs]
SIFbottom = [[K1bottom(c,d) for d in ds] for c in cs]
# create a grid for contouring

DS,CS = np.meshgrid(ds,cs)
"""
fig = plt.figure()
#plt.title('Variable Temperature')
plt.xlabel('Crack Radius ($m$)')
plt.ylabel('Depth Below Surface ($m$)')
plt.ylim(max(ds),min(ds))
plt.xlim(0.01,max(cs))
plt.contourf(CS,DS,SIF,np.linspace(0,np.nanmax(SIF),100),cmap=plt.get_cmap('RdYlBu_r'))
plt.colorbar(ticks=np.arange(0,np.nanmax(SIF),25),label="Stress Intensity Factor $kPa * m^{1/2}$")
plt.hold(True)
plt.contour(CS,DS,SIF,[K1c],colors='k',linewidths=2)
def rho(d):
    return rhoi
plt.contour(CS,DS,[[K1calc(c,d) for d in ds] for c in cs],[K1c],colors='k',linestyles='--',linewidths=2)
plt.savefig('DensityDependent(Integrated).png')#"""

fig = plt.figure(figsize=(9,6))

ax1 = plt.subplot(121)
plt.title('top of crack')
plt.xlabel('Crack Radius ($m$)')
plt.ylabel('Depth Below Surface ($m$)')
plt.ylim(max(ds),min(ds))
plt.xlim(0.01,max(cs))
plt.contourf(CS,DS,SIFtop,np.linspace(np.nanmin(SIFbottom),np.nanmax(SIFbottom),100),cmap=plt.get_cmap('RdYlBu_r'))
#plt.colorbar(ticks=np.arange(0,np.nanmax(SIFtop),25),label="Stress Intensity Factor $kPa * m^{1/2}$")
plt.hold(True)
plt.contour(CS,DS,SIFtop,[K1c],colors='k',linewidths=2)

ax2 = plt.subplot(122,sharey=ax1)
plt.title('bottom of crack')
plt.xlabel('Crack Radius ($m$)')
#plt.ylabel('Depth Below Surface ($m$)')
plt.ylim(max(ds),min(ds))
plt.xlim(0.01,max(cs))
plt.contourf(CS,DS,SIFbottom,np.linspace(np.nanmin(SIFbottom),np.nanmax(SIFbottom),100),cmap=plt.get_cmap('RdYlBu_r'))
plt.colorbar(ticks=np.arange(np.nanmin(SIFbottom),np.nanmax(SIFbottom),10),label="Stress Intensity Factor $kPa * m^{1/2}$")
plt.hold(True)
plt.contour(CS,DS,SIFbottom,[K1c],colors='k',linewidths=2)

plt.savefig('TopBottom.png')