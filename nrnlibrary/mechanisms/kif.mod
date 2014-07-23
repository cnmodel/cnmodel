TITLE kif.mod   DCN pyramidal cell model  
 
COMMENT


This model implements a fast transient potassium current from 
Dorsal Cochlear Nucleus Pyramidal cells
based on kinetic data from Kanold and Manis (1999) and Kanold's dissertation (1999)

-- 15 Jan 1999 P. Manis

2/10/02. P. Manis.

ENDCOMMENT
 
 
UNITS {
 (mA) = (milliamp)
 (mV) = (millivolt)
}
 

NEURON {
 SUFFIX kif
 USEION k READ ek WRITE ik

 
 RANGE gkif, kif_a_inf, kif_i_inf	: fast inactivating potassium current
 RANGE akif, gbar
 RANGE kif_a_tau, kif_i_tau	
 RANGE kif_a_start, kif_i_start
 RANGE  kif_ivh, kif_avh, kif_hivh

 
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
 v (mV)
 celsius (degC)
 dt (ms)
 ek (mV) : = -81.5 (mV)
 ena (mV) : = 50.0 (mV)
 gbar = 0.0125 (mho/cm2) <0,1e9>
 
 kif_ivh = -89.6
 kif_avh = -57.0
 kif_hivh = -87.0
 kif_a_start = -1
 kif_i_start = -1
 }
 
STATE {
        kifa kifi
        }
 
ASSIGNED {
	gkif (mho/cm2)
	ik (mA/cm2)
	kif_a_inf kif_i_inf
	kif_a_tau kif_i_tau
	akif
}

LOCAL kif_a_exp, kif_i_exp
 
? currents
BREAKPOINT {
	SOLVE states METHOD cnexp
	akif = kifa*kifa*kifa*kifa*kifi
	gkif = gbar*akif
	ik = gkif*(v - ek)
}
? currents

UNITSOFF 
 
INITIAL {
	rates(v)
	if(kif_a_start < 0) {		: if xx_(i/a)_start is > 0, then perturbation is done at onset of computations.
		kifa = kif_a_inf
	} else {
		kifa = kif_a_start
	}
	if(kif_i_start < 0) {
		kifi = kif_i_inf
	} else {
		kifi = kif_i_start
	}
}


? states
DERIVATIVE states {  
	rates(v)

	kifa' = (kif_a_inf - kifa) / kif_a_tau
	kifi' = (kif_i_inf - kifi) / kif_i_tau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
	LOCAL  alpha, beta, sum
	TABLE kif_a_inf, kif_a_tau, kif_i_inf, kif_i_tau DEPEND celsius, kif_avh, kif_ivh FROM -200 TO 100 WITH 400

UNITSOFF
	q10 = 3^((celsius - 22)/10)

	: "kif" fast inactivation potassium channel - activation and inactivation 
	kif_a_inf = kif_m(v)
	kif_i_inf = kif_h(v)
	kif_a_tau = kif_mt(v)
	kif_i_tau = kif_ht(v) 

}
 
: Make these as functions so we can view them from hoc, although this 
: may slow things down a bit


FUNCTION kif_m(x) { : ikif activation
	kif_m = 1/(1+exp(-(x-kif_avh)/25.8))
}	

FUNCTION kif_h(x) { : ikif inactivation
	kif_h = 1/(1+exp((x-kif_ivh)/6.7))
}	

FUNCTION kif_mt(x) { : ikif activation tau
	kif_mt = 0.15*exp((x-kif_avh)/10) + 0.3*exp(-(x-kif_avh)/10)
	kif_mt = 0.5 + (1/kif_mt)
}

FUNCTION kif_ht(x) { : ikif inactivation tau
	kif_ht = 0.015*exp((x-kif_hivh)/20)+0.03*exp(-(x-kif_hivh)/20)
	kif_ht = 10 + (1/kif_ht)
}


UNITSON

