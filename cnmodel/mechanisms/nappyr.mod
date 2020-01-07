TITLE nappyr.mod   Persistent sodium conductance
 
COMMENT

Persistent sodium current from Ceballos et al. 2016
Pulled from their file: fuscurr.mod

ENDCOMMENT
 
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}
 

NEURON {
    THREADSAFE
    SUFFIX nappyr
    USEION na READ ena WRITE ina

    RANGE nap_inf, nap_tau
    RANGE gbar, gnap
    RANGE nap

}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
    v (mV)
    celsius (degC)
    dt (ms)
    ena (mV) : = 50.0 (mV)
    gbar = 0.0001 (mho/cm2) <0,1e9>
    q10tau = 3.0
    q10g = 2.0

    nap_vshift = 0 (mV)
    nap_tau = 5.3 (ms)
}
 
STATE {
        nap 
}
 
ASSIGNED {
	gnap (mho/cm2)
	ina (mA/cm2)
	nap_inf
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gnap = gbar*nap*nap*nap
	ina = gnap*(v-ena)
}

 
INITIAL {

	rates(v)
	nap = nap_inf

}


DERIVATIVE states {  
	rates(v)
	nap' = (nap_inf - nap) / nap_tau
}
 

PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.

	: "nap" persistent sodium system
    nap_inf = nap_m(v+nap_vshift)
}
 

FUNCTION nap_m(v (mV)) { : persistent sodium activation
	nap_m = 1.0/(1.0+exp(-(v+55.0)/2.8))
}




