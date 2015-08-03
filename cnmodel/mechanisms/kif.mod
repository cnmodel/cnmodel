TITLE kif.mod   DCN pyramidal cell model fast transient current
    
COMMENT
    
This model implements a fast transient potassium current from 
Dorsal Cochlear Nucleus Pyramidal cells
based on kinetic data from Kanold and Manis (1999) and Kanold's dissertation (1999)
    -- 15 Jan 1999 P. Manis
    2/10/02. P. Manis.
            Pulled from pyr.mod 7/24/2014
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
    RANGE kif_ivh, kif_avh, kif_hivh
}
    
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
    
PARAMETER {
    v (mV)
    celsius (degC)
    dt (ms)
    ek (mV) : = -81.5 (mV)
    ena (mV) : = 50.0 (mV)
    gbar = 0.0125 (mho/cm2) <0,1e9>
    
    kif_ivh = -89.6 (mV)
    kif_avh = -57.0 (mV)
    kif_hivh = -87.0 (mV)
    kif_a_start = -1 (1)
    kif_i_start = -1 (1)
    }
    
STATE {
           kifa kifi
}
    
ASSIGNED {
    gkif (mho/cm2)
    ik (mA/cm2)
    kif_a_inf (1)
    kif_i_inf  (1)
    kif_a_tau (ms)
    kif_i_tau (ms)
    akif ()
    q10 ()
}


BREAKPOINT {
    SOLVE states METHOD cnexp
    akif = kifa*kifa*kifa*kifa*kifi
    gkif = gbar*akif
    ik = gkif*(v - ek)
}

INITIAL {
    q10 = 3^((celsius - 22)/10 (degC))
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
    
DERIVATIVE states {
    rates(v)
    kifa' = (kif_a_inf - kifa) / kif_a_tau
    kifi' = (kif_i_inf - kifi) / kif_i_tau
}
    

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
:    LOCAL  alpha, beta, sum
:    TABLE kif_a_inf, kif_a_tau, kif_i_inf, kif_i_tau DEPEND celsius, kif_avh, kif_ivh FROM -200 TO 100 WITH 400
    : "kif" fast inactivation potassium channel - activation and inactivation
    kif_a_inf = kif_m(v)
    kif_i_inf = kif_h(v)
    kif_a_tau = kif_mt(v)
    kif_i_tau = kif_ht(v)
    }
    
: Make these as functions so we can view them from hoc, although this 
: may slow things down a bit
    
FUNCTION kif_m(v ) { : ikif activation
    kif_m = 1.0/(1+exp(-(v-kif_avh)/25.8 (mV)))
}	

FUNCTION kif_h(v (mV)) { : ikif inactivation
    kif_h = 1.0/(1+exp((v-kif_ivh)/6.7 (mV)))
}

FUNCTION kif_mt(v (mV)) (ms) { : ikif activation tau
    LOCAL x
    x = 0.15 * exp((v-kif_avh)/10 (mV)) + 0.3 *exp(-(v-kif_avh)/10 (mV))
    x = 0.5 + (1.0 /x)
    kif_mt = (x  * 1.0 (ms))/q10
}

FUNCTION kif_ht(v (mV)) (ms) { : ikif inactivation tau
    LOCAL x
    x = 0.015 * exp((v-kif_hivh)/20 (mV))+0.03*exp(-(v-kif_hivh)/20 (mV))
    x = 10 + (1./x)
    kif_ht = (x * 1.0 (ms)) /q10
}


