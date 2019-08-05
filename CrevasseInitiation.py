# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 16:48:54 2015

@author: bh101564

Linear Elastic Fracture Mechanics (LEFM) for varying temperature 
through near surface ice. This shows how it may be possible for
crevasses to initiate beneath the surface. 
"""

import numpy as np
import matplotlib.pyplot as plt

# Thermal Constants (Van der Veen 2013 - from Yen CRREL 1981)
### Constants ###
g = 9.81                            #Gravity (m/s2)
spy  = 60.*60.*24.*365.24           #Seconds per year
gamma = 7.4e-8                      #Clausius-Clapeyron (K Pa-1) From van der Veen p. 209
# pg 144 Van der Veen - from Yen CRREL 1981
rhoi = 917.                         #Bulk density of ice (kg/m3)
kice = 2.1*spy                          #Conductivity of ice (J/mKs)
Cice = 2097.                        #Heat capacity of ice (J/kgK) - ** Van der Veen uses 2097 but see Tr and Aschwanden 2012 and Brinkerhoff 2013
Lf = 3.335e5                        #Latent heat of fusion (J/kg)
# other constants outside of van der Veen
T0 = 273.15                         #Reference Tempearature, triple point for water 
# Strain Enhancement terms
Qc = 7.88e4         # activation energy (J K-1 mol-1) for low trial 6.0e4 for high 11.5e4 for control 7.88e4
A_ = 3.5e-25*spy    # prefactor constant (s-1 Pa-3) for low trial and control 3.5e-25 for high trial 24.e-25
n = 3.0      # creep exponent in Glen's law
R = 8.321           # the gas constant (J K-1 mol-1)
# Fracture Mechanics
K1c = 100.                   # Fracture Toughness (kPa m^.5) Sculson and Duval 2009
edot = .06                 # effective strain rate (yr-1, see Sculson and Duval pg 79; Ryser 2014; Raymond and Malone 1986)
# Geometry-dependent constant, should be of order unity
# Penny Crack (Lawn 1993)
psi = 2./(np.pi)**.5
# from (Tada 1973)
#Y = (1. - 0.025*(a/d)**2 + 0.06*(a/d)**4)*np.sqrt(1/np.cos((np.pi/2.)*(a/d)))

######################################################################################

### Draw a springtime temperature profile based on the equation in Aschwandens 2012 Thermodynmaics course
w = 2*np.pi                 # period of surface oscillation
ds = np.linspace(0,10,100)  # depth below the surface (m)
MoY = 2.                    # Month of year for the desired temperature profile
Tc = T0-10.0                  # Temperature profile oscillates around the convergence temp T0 (K)
Ta = 10.0                   # Size of annual oscillations at the surface (K)

def VarTemp(d):
    return Tc + Ta*np.exp(-d*np.sqrt(w/(2.*kice/(rhoi*Cice))))*np.sin((w*(MoY/12.)) - d*np.sqrt(w/(2.*kice/(rhoi*Cice))))

### Calculate the rate factor , A, over a range of temperatures based on the equation in Van der Veen 2013 (p 33)
def RateFactor(T):
    return A_*np.exp(-(Qc/R)*((1./(T))-(1./(T0-10.))))  # rate factor (Pa-3 yr-1)  

"""
### Calculate the viscosity parameter, B, over a range of temperatures based on the equation in Van der Veen 2013 (p 33)
B0 = 2.207e-3               # Viscosity constant (kPa yr^1/3)                     
T0 = 3155                   # Temperature constant (K)

def ViscosityPar(T):
    return B0*np.exp((T0/T)-(C/((Tr-T)**K))) #Viscosity Parameter (kPa yr^1/3)
    
"""

######################################################################################
### Fracture Mechanics  

fig = plt.figure()
ax1 = plt.subplot(121)
plt.plot
plt.ylabel('Depth ($m$)')
plt.xlabel('Temperature ($^\circ C$)')
plt.ylim(10,0)
plt.plot(VarTemp(ds)-T0,ds,'k',lw=2)
plt.plot(np.ones_like(ds)*Tc-T0,ds,'k--',lw=2)
plt.grid()

ax2 = plt.subplot(122)
plt.tick_params(axis='both',which='both',labelleft='off')

c = 0.2
sigma = ((edot/RateFactor(VarTemp(ds)))**(1/n) - rhoi*g*ds)/1000. 
# Stress Intensity Factor (kPa m**1/2)
VarK1s = ((edot/RateFactor(VarTemp(ds)))**(1/n) - rhoi*g*ds)/1000. *np.sqrt(np.pi*c)
K1s = ((edot/RateFactor(Tc))**(1/n) - rhoi*g*ds)/1000. *np.sqrt(np.pi*c)

plt.fill_betweenx(np.linspace(-10,10,2),K1c,200,color='k',alpha=0.15)
plt.plot(VarK1s,ds,'k',lw=2)
plt.plot(K1s,ds,'k--',lw=2)
plt.grid()
plt.ylim(10,0)
plt.xlim(round(min(K1s),-1),round(max(K1s),-1))
plt.xlabel('Stress Intensity Factor ($kPa\ m^{1/2}$)')
plt.savefig('StressIntensityFactors.png')

"""
# Calculate the effective stress based on some constant effective strain rate

Rxx = (edot/RateFactor(VarTemp(ds)))**(1/n)
Lxx = -rhoi*g*ds 

# Find the Stress Intensity factor based on the total applied stress
cs = np.linspace(0,.4,100)       # crack sizes = 2a (m)
K1s = np.zeros((len(cs),len(ds)))                    # Stress Intensity Factor
VarK1s = np.zeros((len(cs),len(ds)))                    # Stress Intensity Factor (variable temps)
i=0
for c in cs:
    j=0
    for d in ds:
        # Total applied stress (kPa)
        sigma = ((edot/RateFactor(VarTemp(d)))**(1/n) - rhoi*g*d)/1000. 
        # Stress Intensity Factor (kPa m**1/2)
        VarK1s[i,j] = sigma*np.sqrt(np.pi*c)
        K1s[i,j] = ((edot/RateFactor(Tc))**(1/n) - rhoi*g*d)/1000. *np.sqrt(np.pi*c)
        j+=1    
    i+=1

D,C = np.meshgrid(ds,cs)

fig = plt.figure()
ax1 = plt.subplot(121)
plt.plot
plt.ylabel('Depth ($m$)')
plt.xlabel('Temperature ($^\circ C$)')
plt.ylim(10,0)
plt.plot(VarTemp(ds)-T0,ds,'k',lw=2)
plt.plot(np.ones_like(ds)*Tc-T0,ds,'k--',lw=2)
plt.grid()

ax2 = plt.subplot(122)
plt.tick_params(axis='both',which='both',labelleft='off')
plt.xlabel('Starter Crack Length ($m$)')
ax2.set_xticks(np.arange(0,0.5,0.1))
plt.ylim(10,0)
#plt.contourf(C,D,VarK1s,np.linspace(0,200,100),cmap=plt.get_cmap('RdYlBu_r'))
#plt.colorbar(ticks=np.arange(0,201,25),label="Stress Intensity Factor ($kPa\ m^{1/2})$")
#plt.hold(True)
plt.contour(C,D,VarK1s,[K1c],colors='k',linewidths=2)
plt.contour(C,D,K1s,[K1c],colors='k',linestyles='--',linewidths=2)
plt.grid()
plt.savefig('StressIntensityFactors.png')#"""