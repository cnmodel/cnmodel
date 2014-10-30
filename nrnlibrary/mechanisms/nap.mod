TITLE nap.mod   Persistent sodium conductance
 
COMMENT


Persistent sodium current from deSchutter and Bower, J. Neurophys.
71:375, 1994.


2/10/02, 10/10/2014. P. Manis.

ENDCOMMENT
 
 
UNITS {
 (mA) = (milliamp)
 (mV) = (millivolt)
}
 

NEURON {
 THREADSAFE
 SUFFIX nap
 USEION na READ ena WRITE ina
 
 RANGE nap_inf, nap_tau, napi_inf, napi_tau
 RANGE gbar, gnap
 RANGE nap_A, nap_B, nap_C, nap_D, nap_E, nap_F, nap_G, nap_H
 
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
 v (mV)
 celsius (degC)
 dt (ms)
 ena (mV) : = 50.0 (mV)
 gbar = 0.00001 (mho/cm2) <0,1e9>
 q10tau = 3.0
 q10g = 2.0
 
 nap_shift = 0 (mV)
 
 nap_A = 200  (/ms): parameters from Bowers et al.
 nap_B = 1
 nap_C = -18 (mV)
 nap_D = -16 (mV)
 nap_E = 25 (/ms)
 nap_F = 1
 nap_G = 58 (mV)
 nap_H = 8 (mV)
 }
 
STATE {
        nap napi 
        }
 
ASSIGNED {
	gnap (mho/cm2)
	ina (mA/cm2)
	nap_inf
	nap_tau (ms)
    nap_tau1 (/ms)
    nap_tau2 (/ms)
	napi_inf
	napi_tau (ms)
    napi_tau1 (/ms)
    napi_tau2 (/ms)
    qg ()  : computed q10 for gnabar based on q10g
    q10 ()
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gnap = gbar*nap*nap*nap*napi
	ina = gnap*(v-ena)
}

 
INITIAL {
    qg = q10g^((celsius-22)/10 (degC))
    q10 = q10tau^((celsius - 22)/10 (degC)) : if you don't like room temp, it can be changed!
	rates(v)
	nap = nap_inf
	napi = napi_inf
}


DERIVATIVE states {  
	rates(v)
	nap' = (nap_inf - nap) / nap_tau
	napi' = (napi_inf - napi) / napi_tau
}
 

PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
    LOCAL x
	: "nap" persistent sodium system
	nap_inf = na_p(v + nap_shift)
	nap_tau1 = nap_A/(nap_B + exp((v + nap_shift + nap_C)/nap_D))
    nap_tau2 = nap_E/(nap_F + exp((v + nap_shift + nap_G)/nap_H))
	nap_tau = 1./(nap_tau1 + nap_tau2)

	:nap_tau = na_ptau(v + nap_shift)

	: "nap" persistent sodium system - inactivation...
	napi_inf = na_pi(v + nap_shift)
    napi_tau1 = (1 (/ms)) /(0.06435/(1+exp((v + nap_shift + 73.26415)/3.71928 (mV))))
    napi_tau2 = (0.13496 (/ms))/(1 +exp((v + nap_shift + 10.27853)/(-9.09334 (mV))))
	napi_tau = 1 /(napi_tau1 + napi_tau2)
}
 
LOCAL p

FUNCTION na_p(v (mV)) { : persistent sodium activation
: Bowers
	p = nap_A/(nap_B + exp((v + nap_shift + nap_C)/nap_D))
	na_p = p/(p+nap_E/(nap_F + exp((v + nap_shift + nap_G)/nap_H)))
}


FUNCTION na_pi(x (mV)) { : persistent sodium inactivation
: Bowers
	na_pi = 0.06435/(1+exp((x+73.26415)/3.71928 (mV)))
	na_pi = na_pi/(na_pi + (0.13496/(1+exp((v+10.27853)/(-9.09334 (mV))))))
}


