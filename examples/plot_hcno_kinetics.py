#!/usr/bin/python
"""
test and plot hcno stuff

"""
import numpy as np
import matplotlib.pyplot as plt

#Parameters#
gbar = 0.0005       # (mho/cm2)   
                            
vhalf1  = -50   # (mV)        # v 1/2 for forward
vhalf2  = -84   # (mV)        # v 1/2 for backward    
gm1   = 0.3 ## (mV)           # slope for forward
gm2   = 0.6  #    # (mV)      # slope for backward
zeta1   = 3 #   (/ms)       
#    zeta2   = 3 #   (/ms)       
a01 = 0.008  #(/ms)
a02 = 0.0029 #(/ms)
frac = 0.0
c0 = 273.16  #(degC)
thinf  = -66    # (mV)        # inact inf slope   
qinf  = 7   # (mV)        # inact inf slope 
q10tau = 4.5             # from Magee (1998)
#v       # (mV)
q10g = 4.5    # Cao and Oertel
celsius = 35.

F = 9.648e4
R = 8.314

ct = 1e-3*zeta1*F/(R*(c0+celsius))

q10 = q10tau**((celsius - 33.0)/10.0) # (degC)) : if you don't like room temp, it can be changed!

def rates(v): # (mV)) {  
    tau1 = bet1(v)/(q10*a01*(1.0+alp1(v)))
    tau2 = bet2(v)/(q10*a02*(1.0+alp2(v)))
    hinf = 1.0/(1.0+np.exp((v-thinf)/qinf))
    return tau1, tau2, hinf

def alp1(v):# (mV)) {
    alp1 = np.exp((v-vhalf1)*ct) 
    return alp1

def bet1(v):# (mV)) {
    bet1 = np.exp(gm1*(v-vhalf1)*ct)
    return bet1

def alp2(v):# (mV)) {
    alp2 = np.exp((v-vhalf2)*ct) 
    return alp2

def bet2(v):# (mV)) {
    bet2 = np.exp(gm2*(v-vhalf2)*ct) 
    return bet2

v = np.linspace(-120., 0., 120)

t1, t2, hi = rates(v)

print v
plt.figure(1)
plt.plot(v, t1)
plt.plot(v, t2)
plt.figure(2)
plt.plot(v, hi)
plt.show()
