TITLE ihpyr.mod   DCN pyramidal cell model H-current
 
COMMENT

This model is part a Dorsal Cochlear Nucleus Pyramidal point cell
based on kinetic data from Kanold and Manis (1999) and Kanold's dissertation (1999)

-- 15 Jan 1999 P. Manis

Added export of start states for some variables to do perturbation tests
These start values replace the "inf" values used in the initialization procedure
Note that if the start variable is set to a value less than 0, 
then the default initialization will be done. Typically I use a value of -1 for this flagging
Note also that it is possible to set the initial values > 1 but this is meaningless in terms of
the present equations. 
-- 5 Feb 1999 P. Manis

Added Patrick's version of ih as ihpyr
Model is from Destexhe and Babloyantz 1993; Destexhe et al. 1993


2/10/02. P. Manis.
7/23/2014 P. Manis - separated from pyr.mod.

ENDCOMMENT
 
 
UNITS {
 (mA) = (milliamp)
 (mV) = (millivolt)
}
 

NEURON {
 THREADSAFE
 SUFFIX ihpyr
 USEION na READ ena WRITE ina
 USEION k READ ek WRITE ik
 NONSPECIFIC_CURRENT i
: USEION h READ eh WRITE ih VALENCE 1
 RANGE eh
:
 RANGE gh, kh_m_inf, kh_n_inf, aih, gbar, ghvshift
 RANGE kh_m_tau, kh_n_tau

}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
 v (mV)
 celsius (degC)
 dt (ms)
 ek (mV) : = -81.5 (mV)
 ena (mV) : = 50.0 (mV)
 gbar = 0.00025 (mho/cm2) <0,1e9>
 ghvshift = 0 (mV)
 eh (mV) : = -43.0(mV)
}
 
STATE {
       khm khn
}
 
ASSIGNED {
	gh (mho/cm2)
	ina (mA/cm2)
	ik (mA/cm2)
	ih (mA/cm2)
    i (mA/cm2)
	kh_m_inf kh_n_inf
	kh_m_tau kh_n_tau	
	aih
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	aih = khm*khn
	gh = gbar*aih
 	ih = gh*(v - eh)
    i = ih
}

UNITSOFF 
 
INITIAL {
	rates(v)
	khm = kh_m_inf
	khn = kh_n_inf
}


DERIVATIVE states {  
	rates(v)
	khm' = (kh_m_inf - khm) / kh_m_tau
	khn' = (kh_n_inf - khn) / kh_n_tau
}


LOCAL q10

PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
	                      :Call once from HOC to initialize inf at resting v.
	q10 = 3^((celsius - 22)/10 (degC))

	:"kh" adaptation of Destexhe hyp-activated cation current by Patrick Kanold
	kh_m_inf = kh_m(v) 
	kh_n_inf = kh_n(v)
	kh_m_tau = kh_mt(v)
	kh_n_tau = kh_nt(v) 
}
 
: Make these as functions so we can view them from hoc, although this 
: may slow things down a bit


FUNCTION kh_m(x (mV)) {
	kh_m = 1/(1+exp((x+68.9+ghvshift)/6.5))
}

FUNCTION kh_n(x (mV)) {
	kh_n = 1/(1+exp((x+68.9+ghvshift)/6.5)) : same as kh_m, but for completeness, compute this
}

FUNCTION kh_mt(v (mV)) {
	kh_mt = exp((v+183.6+ghvshift)/15.24)
}

FUNCTION kh_nt(v (mV)) {
	kh_nt = exp((v+158.6+ghvshift)/11.2)/(1+exp((v+75+ghvshift)/5.5))
}



