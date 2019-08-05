# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 10:57:30 2016

@author: ben
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

### Stress Intensity Factor for a constant stress
ds = np.arange(100.)
H1,H2 = 2000.,500.          # Ice sheet thickness (m)
Rxx = 150.                   # Dynamic Tensile Stress (kPa)
def F(d,H):                 # Geometric constant
    l = d/H    
    return 1.12 - 0.23*l + 10.55*l**2 - 21.72*l**3 + 30.39*l**4
K1 = np.array([[F(d,H)*Rxx*np.sqrt(np.pi*d)/1000. for d in ds] for H in [H1,H2]])

"""
fig = plt.figure()
plt.title('Stress Intensity Factor for Rxx=150kPa')
plt.xlim(0,4)
plt.ylim(100,0)
plt.xlabel("$K_1 (MPa * m^{1/2})$")
plt.ylabel("d($m$)")
plt.plot(K1[0],ds,'k-',lw=3,label='2000 $m$')
plt.plot(K1[1],ds,'k-',lw=2,label='500 $m$')
plt.legend()
plt.savefig('Rxx_StressIntensityFactor.png')"""


### Stress Intensity Factor for a depth varying stress
g = 9.81                    #Gravity (m/s2)
rhoi = 917.                 #Bulk density of ice (kg/m3)
rhos = 350.         # Density of surface snow
C = 0.02            # Densification rate
H = 2000.
def L(b):
    return -rhoi*g*b + ((rhoi-rhos)/C)*g*(1-np.exp(-C*b))
def G(b,d,H):
    gam = b/d
    lam = d/H
    return (3.52*(1-gam))/((1-lam)**1.5) - (4.35-5.28*gam)/((1-lam)**.5) +\
            ((1.3 - 0.3*gam**1.5)/((1-gam**2)**.5) + 0.83 - 1.76*gam)*(1-(1-gam)*lam)
def Func(b):
    return L(b)*G(b,d,H)
K2 = [[2./(np.sqrt(np.pi*d))*(quad(Func,0,d)[0]) for d in ds]]

def L(b):
    return -rhoi*g*b
K2 = np.append(K2,[[2./(np.sqrt(np.pi*d))*(quad(Func,0,d)[0]) for d in ds]],axis=0)
K2 = K2/1000000.


fig = plt.figure()
plt.plot()
plt.title('Stress Intensity Factor from Lithostatic Stress 2000m Glacier')
plt.xlim(-12,0)
plt.ylim(100,0)
plt.xlabel("$K_2 (MPa * m^{1/2})$")
plt.ylabel("d($m$)")
plt.plot(K2[1],ds,'k-',lw=2,label='Constant Density')
plt.plot(K2[0],ds,'k-',lw=3,label='Variable')
plt.legend(loc = 2)
#plt.savefig('L_StressIntensityFactor.png')
"""


### Add together the Lithostatic and the Tensile stress intensity factors for two different stress regimes
K1 = np.array([[F(d,H1)*R*np.sqrt(np.pi*d)/1000. for d in ds] for R in [100,200]])

K1tot = [K+K2[0] for K in K1]

fig = plt.figure()
plt.plot()
plt.title('Net Stress Intensity Factor')
plt.xlim(0,1.6)
plt.ylim(60,0)
plt.xlabel("$K_net (MPa * m^{1/2})$")
plt.ylabel("d($m$)")
plt.axvline(0.1,linestyle='--',color='k')
plt.axvline(0.4,linestyle='--',color='k')
plt.plot(K1tot[0],ds,'k-',lw=2,label='100 kPa')
plt.plot(K1tot[1],ds,'k-',lw=3,label='200 kPa')
plt.legend()
plt.savefig('Total_StressIntensityFactor.png')#"""