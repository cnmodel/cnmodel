TITLE pyr.mod   DCN pyramidal cell model  
 
COMMENT

Revised version of DCN Pyramidal cell model based on new hh.hoc file in NEURON

This model implements a Dorsal Cochlear Nucleus Pyramidal point cell
based on kinetic data from Kanold and Manis (1999) and Kanold's dissertation (1999)

-- 15 Jan 1999 P. Manis

Added export of start states for some variables to do perturbation tests
These start values replace the "inf" values used in the initialization procedure
Note that if the start variable is set to a value less than 0, 
then the default initialization will be done. Typically I use a value of -1 for this flagging
Note also that it is possible to set the initial values > 1 but this is meaningless in terms of
the present equations. 
-- 5 Feb 1999 P. Manis

Removed slow Ih current 30 Jan 2000. P. Manis
- also renamed variables to saner forms
Added Patrick's version of ih as ihd

Added persistent sodium current from deSchutter and Bower, J. Neurophys.
71:375, 1994.
Changed to Nap current from Quadroni and Knopfel, 1994


2/10/02. P. Manis.

ENDCOMMENT
 
 
UNITS {
 (mA) = (milliamp)
 (mV) = (millivolt)
}
 

NEURON {
 SUFFIX pyr
 USEION na READ ena WRITE ina
 USEION k READ ek WRITE ik
 USEION h READ eh WRITE ih VALENCE 1
 NONSPECIFIC_CURRENT il
 
 RANGE gna, gk, minf, hinf, ninf, gnabar, gkbar		: sodium channels and delayed rectifier
 RANGE mtau, htau, ntau: time constants for sodium channels and delayed rectifier
 RANGE kd_avh
 
 RANGE nap_inf, nap_tau, napi_inf, napi_tau
 RANGE gnapbar
 RANGE nap_A, nap_B, nap_C, nap_D, nap_E, nap_F, nap_G, nap_H
 
 RANGE gkif, kif_a_inf, kif_i_inf	: fast inactivating potassium current
 RANGE akif, gkifbar
 RANGE kif_a_tau, kif_i_tau	
 RANGE kif_a_start, kif_i_start
 
 RANGE gkis, kis_a_inf, kis_i_inf, gkisbar	: slow inactivating potassium current
 RANGE kis_a_tau, kis_i_tau	
 RANGE kis_a_start, kis_i_start 	:flags for perturbation analysis
 RANGE akis				: state variable, slow
 RANGE kif_avh, kif_ivh, kis_avh, kis_ivh : adjustable ikif and
			: ikis voltage dep...

 RANGE gh, kh_m_inf, kh_n_inf, aih, ghbar, ghvshift
 RANGE kh_m_tau, kh_n_tau

 RANGE gl, el							: leak (linear, ohmic with reversal at el)

}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
 v (mV)
 celsius (degC)
 dt (ms)
 ek (mV) : = -81.5 (mV)
 ena (mV) : = 50.0 (mV)
 gnabar =  0.02857 (mho/cm2)	<0,1e9>
 gnapbar = 0.00001 (mho/cm2) <0,1e9>
 gkbar = 0.006667 (mho/cm2)	<0,1e9>
 gkifbar = 0.0125 (mho/cm2) <0,1e9>
 gkisbar = 0.0033333 (mho/cm2)
 ghbar = 0.00025 (mho/cm2) <0,1e9>
 ghvshift = 0 (mV)
 eh (mV) : = -43.0(mV)
 gl = 0.000250 (mho/cm2)	<0,1e9>
 el = -57.7 (mV) 
 mtau0 = 0.05 (ms) <0.01,100>
 htau0 = 0.5 (ms) <0.1,100>
 ntau = 0.5 (ms) <0.1,100>
 kd_avh = -40
 nap_shift = 0
 
 nap_A = 200  : parameters from Bowers et al.
 nap_B = 1
 nap_C = -18
 nap_D = -16
 nap_E = 25
 nap_F = 1
 nap_G = 58
 nap_H = 8
 
 
 kif_ivh = -89.6
 kif_avh = -57.0
 kif_hivh = -87.0
 kis_ivh = -40.9
 kis_avh = -38.4
 kif_a_start = -1
 kif_i_start = -1
 kis_a_start = -1
 kis_i_start = -1
 }
 
STATE {
        m h n kifa kifi kisa kisi khm khn nap napi 
        }
 
ASSIGNED {
	gna (mho/cm2)
:	gnap (mho/cm2) : bowers
	gnap (mho/cm2) : knopfel
	gk (mho/cm2)
	gkif (mho/cm2)
	gkis (mho/cm2)
	gh (mho/cm2)
	ina (mA/cm2)
	ik (mA/cm2)
	ih (mA/cm2)
	il (mA/cm2)
	minf hinf mtau htau
	ninf 
	nap_inf nap_tau
	napi_inf napi_tau
	kif_a_inf kif_i_inf
	kis_a_inf kis_i_inf
	kif_a_tau kif_i_tau
	kis_a_tau kis_i_tau
	kh_m_inf kh_n_inf
	kh_m_tau kh_n_tau	
	akif
	akis
	aih
}

LOCAL mexp, hexp, nexp, nap_exp, napi_exp, kif_a_exp, kif_i_exp, kis_a_exp, kis_i_exp, kh_m_exp, kh_n_exp
 
? currents
BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gnabar*m*m*h
	gnap = gnapbar*nap*nap*nap*napi
:	gnap = gnapbar * nap
	ina = gna*(v - ena) + gnap*(v-ena) + 0.25*gh*(v-ena)
	gk = gkbar*n*n
	akif = kifa*kifa*kifa*kifa*kifi
	gkif = gkifbar*akif
	akis = kisa*kisa*kisa*kisa*kisi
	gkis = gkisbar*akis
	aih = khm*khn
	gh = ghbar*aih
	ik = gk*(v - ek) + gkif*(v - ek) + gkis*(v - ek) + 0.75*gh*(v-ek)
 	ih = gh*(v - eh)
	il = gl*(v - el)
}
? currents

UNITSOFF 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
	nap = nap_inf
	napi = napi_inf
	n = ninf
	khm = kh_m_inf
	khn = kh_n_inf
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
	if(kis_a_start < 0) {
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


? states
DERIVATIVE states {  
	rates(v)
	m' = (minf - m) / mtau
	h' = (hinf - h) / htau
	nap' = (nap_inf - nap) / nap_tau
	napi' = (napi_inf - napi) / napi_tau
	n' = (ninf - n) / ntau
	kifa' = (kif_a_inf - kifa) / kif_a_tau
	kifi' = (kif_i_inf - kifi) / kif_i_tau
	kisa' = (kis_a_inf - kisa) / kis_a_tau
	kisi' = (kis_i_inf - kisi) / kis_i_tau
	khm' = (kh_m_inf - khm) / kh_m_tau
	khn' = (kh_n_inf - khn) / kh_n_tau
}
 
LOCAL q10


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
	LOCAL  alpha, beta, sum
	TABLE minf, mtau, hinf, htau, nap_inf, nap_tau, napi_inf, napi_tau, ninf, ntau, kif_a_inf, kif_a_tau, kif_i_inf, kif_i_tau, kis_a_inf, kis_a_tau, kis_i_inf, kis_i_tau, kh_m_inf, kh_n_inf, kh_m_tau, kh_n_tau DEPEND celsius, kif_avh, kif_ivh, kis_avh, kis_ivh, kd_avh, nap_A, nap_C, nap_D, nap_E, nap_G, nap_H FROM -200 TO 100 WITH 400

UNITSOFF
	q10 = 3^((celsius - 22)/10)

	: "m" sodium activation system
	minf = na_m(v)
	mtau = na_mt(v)
	
	: "h" sodium inactivation system
	hinf = na_h(v)
	htau = na_ht(v)
	
	: "nap" persistent sodium system
	nap_inf = na_p(v+nap_shift)
	nap_tau = na_pt(v+nap_shift)
	
	: "nap" persistent sodium system - inactivation...
	napi_inf = na_pi(v+nap_shift)
	napi_tau = na_pit(v+nap_shift)
	
	: "n" potassium activation system
    ninf = kd_m(v)

	: "kif" fast inactivation potassium channel - activation and inactivation 
	kif_a_inf = kif_m(v)
	kif_i_inf = kif_h(v)
	kif_a_tau = kif_mt(v)
	kif_i_tau = kif_ht(v) 

	:"kis" slow inactivating potassium channel - activation and inactivation
	kis_a_inf = kis_m(v)
	kis_i_inf = kis_h(v)
	kis_a_tau = kis_mt(v)
	kis_i_tau = kis_ht(v)

	:"kh" adaptation of Destexhe hyp-activated cation current by Patrick Kanold
	kh_m_inf = kh_m(v) 
	kh_n_inf = kh_n(v)
	kh_m_tau = kh_mt(v)
	kh_n_tau = kh_nt(v) 
}
 
: Make these as functions so we can view them from hoc, although this 
: may slow things down a bit

FUNCTION na_m(x) { : sodium activation
	na_m = 1/(1+exp(-(x+38)/3.0))	: POK version
:	na_m = alphbet(x,35,0,5,-10) :de Schutter (doesn't work well in our version)
:	na_m = na_m/(na_m + alphbet(x,7,0,65,20))
}

FUNCTION na_mt(x) { : sodium activation with taus
	na_mt = mtau0 : flat time constants
:	na_mt = alphbet(x,35,0,5,-10)
:	na_mt = 1/(na_mt + alphbet(x,7,0,65,20))
}

FUNCTION na_h(x) { : sodium inactivation
	na_h = 1/(1+exp((x+43)/3.0))	: flat time constants (POK version)
:	na_h = alphbet(x,0.225,1,80,10)	
:	na_h = na_h/(na_h + alphbet(x,7.5,0,-3,-18))
}

FUNCTION na_ht(x) { : sodium inactivation tau
	na_ht = htau0 : POK: flat time constants
:	na_ht = alphbet(x,0.225,1,80,10) : de Schutter version (doesn't work well with other stuff)
:	na_ht = 1/(na_ht + alphbet(x,7.5,0,-3,-18))
}

FUNCTION na_p(x) { : persistent sodium activation
: BOwers
	na_p = alphbet(x,nap_A,nap_B,nap_C,nap_D)
	na_p = na_p/(na_p+alphbet(x,nap_E,nap_F,nap_G,nap_H))
: Knopfel
:	na_p = (0.12*exp(0.12*(x+56)))/( (0.12*exp(0.12*(x+56))) + (0.12*exp(-0.03*(x + 56))))
: Original
:	na_p = 1.032/(1+exp((x+6.99)/-14.87115))
:	na_p = na_p/(na_p + (5.79/(1+exp((v+130.4)/22.9))))
}

FUNCTION na_pt(x) { :persistent sodium tau.
: Bowers
	na_pt = alphbet(x,nap_A,nap_B,nap_C,nap_D)
	na_pt = 1/(na_pt+alphbet(x,nap_E,nap_F,nap_G,nap_H))
: Knopfel 1/alpha+beta 
:	na_pt = (0.12*exp(0.12*(x+56))) + (0.12*exp(-0.03*(x + 56)))
:	na_pt = 1/na_pt 
:	if(na_pt < 1) { : clamp minimum activation rate to 1 msec
:		na_pt = 1
:	}
	: original
:	na_p = 1.032/(1+exp((x+6.99)/-14.87115))
:	na_p = na_p/(na_p + (5.79/(1+exp((v+130.4)/22.9))))

}

FUNCTION na_pi(x) { : persistent sodium inactivation
: Bowers
	na_pi = 0.06435/(1+exp((x+73.26415)/3.71928))
	na_pi = na_pi/(na_pi + (0.13496/(1+exp((v+10.27853)/(-9.09334)))))
: Knopfel
:	na_pi = 1 : not used
}

FUNCTION na_pit(x) { :persistent sodium tau inactivation.
: Bowers
	na_pit = 0.06435/(1+exp((x+73.26415)/3.71928))
	na_pit = 1/(na_pit + (0.13496/(1+exp((v+10.27853)/(-9.09334)))))
:Knopfel
:	na_pit = 1 : not used
}

FUNCTION kd_m(x) { : potassium activation
	kd_m = 1/(1+exp(-(x-kd_avh)/3))		: flat time constants
}

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

FUNCTION kis_m(x) { : ikis activation
	kis_m = 1/(1+exp(-(x-kis_avh)/23.7))
}

FUNCTION kis_h(x) { : ikis inactivation
	kis_h = 1/(1+exp((x-kis_ivh)/9))
}

FUNCTION kis_mt(x) { : ikis activation tau
	kis_mt = 0.15*exp((x-kis_avh)/10) + 0.3*exp(-(x-kis_avh)/10)
	kis_mt = 0.5 + 1/kis_mt
}

FUNCTION kis_ht(x) { : ikis inactivation tau
	kis_ht = 200
}

FUNCTION kh_m(x) {
	kh_m = 1/(1+exp((x+68.9+ghvshift)/6.5))
}

FUNCTION kh_n(x) {
	kh_n = 1/(1+exp((x+68.9+ghvshift)/6.5)) : same as kh_m, but for completeness, compute this
}

FUNCTION kh_mt(v) {
	kh_mt = exp((v+183.6+ghvshift)/15.24)
}

FUNCTION kh_nt(v) {
	kh_nt = exp((v+158.6+ghvshift)/11.2)/(1+exp((v+75+ghvshift)/5.5))
}


FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

FUNCTION alphbet(x,A,B,C,D) { 	: alpha/beta general functions for
				: transcrbing GENESIS models
alphbet = A/(B+exp((x+C)/D))
}

UNITSON

