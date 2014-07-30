TITLE kis.mod   DCN pyramidal cell model Slow transient K current

COMMENT


This model implements the slow transient potassium current from
Dorsal Cochlear Nucleus Pyramidal cells
based on kinetic data from Kanold and Manis (1999) and Kanold's dissertation (1999)

-- 15 Jan 1999 P. Manis

2/10/02, 7/24/2014. P. Manis.

ENDCOMMENT


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}


NEURON {
    THREADSAFE
    SUFFIX kis
    USEION k READ ek WRITE ik
    RANGE gkis, kis_a_inf, kis_i_inf	: fast inactivating potassium current
    RANGE akis, gbar
    RANGE kis_a_tau, kis_i_tau	
    RANGE kis_a_start, kis_i_start
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
    v (mV)
    celsius (degC)
    dt (ms)
    ek (mV) : = -81.5 (mV)
    ena (mV) : = 50.0 (mV)
    gbar = 0.0033333  (mho/cm2) <0,1e9>

    kis_ivh = -40.9 (mV)
    kis_avh = -38.4 (mV)
    kis_a_start = -1
    kis_i_start = -1
}

STATE {
    kisa kisi
}

ASSIGNED {
	gkis (mho/cm2)
	ik (mA/cm2)
	kis_a_inf kis_i_inf
	kis_a_tau (ms)
    kis_i_tau (ms)
	akis
    q10 ()
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	akis = kisa*kisa*kisa*kisa*kisi
	gkis = gbar*akis
	ik = gkis*(v - ek)
}

INITIAL {
	q10 = 3^((celsius - 22)/10 (degC))
	rates(v)
	if(kis_a_start < 0) {		: if xx_(i/a)_start is > 0, then perturbation is done at onset of computations.
		kisa = kis_a_inf
	} else {
		kisa = kis_a_start
	}
	if(kis_i_start < 0) {
		kisi = kis_i_inf
	} else {
		kisi = kis_i_start
	}
}


DERIVATIVE states {
	rates(v)
	kisa' = (kis_a_inf - kisa) / kis_a_tau
	kisi' = (kis_i_inf - kisi) / kis_i_tau
}

PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.

	: "kis" fast inactivation potassium channel - activation and inactivation
	kis_a_inf = kis_m(v)
	kis_i_inf = kis_h(v)
	kis_a_tau = kis_mt(v)
	kis_i_tau = kis_ht(v)
}

: Make these as functions so we can view them from hoc, although this
: may slow things down a bit

FUNCTION kis_m(v (mV)) { : ikis activation
	kis_m = 1/(1+exp(-(v-kis_avh)/23.7 (mV)))
}

FUNCTION kis_h(v (mV)) { : ikis inactivation
	kis_h = 1/(1+exp((v-kis_ivh)/9 (mV)))
}

FUNCTION kis_mt(v (mV)) (ms) { : ikis activation tau
	LOCAL x
    x = 0.15*exp((v-kis_avh)/10 (mV)) + 0.3*exp(-(v-kis_avh)/10 (mV))
    x = 0.5 + (1.0 /x)
    kis_mt = (x  * 1.0 (ms))/q10
}

FUNCTION kis_ht(v (mV)) (ms) { : ikis inactivation tau
	kis_ht = 200 (ms)
}

