TITLE hasenstaub.mod A HH-style model for exploring Na/K fluxes

COMMENT

NEURON implementation of the model from Hasenstaub et al., PNAS, 2010

Contact: pmanis@med.unc.edu

ENDCOMMENT

UNITS {
        (molar) = (/liter)
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
        (microS) = (micromho)
        (mM) = (millimolar)
}

NEURON {
        SUFFIX hstb
        USEION na READ nai, nao, ena WRITE ina
        USEION k READ ki, ko, ek WRITE ik
        NONSPECIFIC_CURRENT il
        RANGE gnabar, gkbar, gl, el, il

}

PARAMETER {
	v (mV)
	celsius (degC) : 20
	gnabar= 7.0 (mho/cm2)
	gkbar = 7.0 (mho/cm2)
	nai (mM) : 6
	nao (mM) : 140
	ki (mM) : 140
	ko (mM) : 2.5
	gl = 0.1 (mho/cm2)
	el = -70 (mV)
	ena = 50 (mV)
	ek = -90 (mV)
}

STATE {
        m h n
}

ASSIGNED {
    ik (mA/cm2)
    ina (mA/cm2)
    il (mA/cm2)
    ninf (1)
    minf (1)
    hinf (1)
    ntau (ms)
    htau (ms)
    mtau (ms)
}


BREAKPOINT {
    LOCAL gna, gk
	SOLVE states METHOD cnexp

	gna = gnabar*(m*m*m)*h
    ina = gnabar*(v - ena)
    gk = gkbar*n*n
    ik = gkbar * (v - ek)
    il = gl * (v - el)
}


INITIAL {
    trates(v)
    m = minf
    h = hinf
    n = ninf
}

DERIVATIVE states {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	m' = (minf - m)/mtau
	h' = (hinf - h)/htau
	n' = (ninf - n)/ntau
}

UNITSOFF

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
    LOCAL am, bm, ah, bh, an, bn
    
    am = -40.0 (ms) *(-75.5+v)/(exp((-75.5+v)/-13.5)-1)
    bm = 1.2262 (ms) /(exp(v/42.24))
    mtau = 1./(am+bm)
    minf = am*mtau
    
    ah = 0.00315 (ms) /(exp(v/24.186))
    bh = -0.0153 (ms) *(51.2+v)/(exp(-(51.2+v)/5.2)-1)
    htau = 1./(ah+bh)
    hinf = ah*htau

    an = -1 (ms) * (-95+v)/(exp((-95+v)/-11.8)-1)
    bn = 0.025/(exp(v/22.222))
    ntau = 1./(an+bn)
    ninf = an*ntau
}

PROCEDURE trates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	TABLE minf, hinf, ninf
	DEPEND celsius FROM -150 TO 150 WITH 300

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc
	}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
UNITSON

