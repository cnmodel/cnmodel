TITLE jsrnaf.mod VCN Na conductance, fast model

COMMENT
gnaf is the modified form used in his
1993 M.S. thesis (as in Rothman et al., J. Neurophysiol. 70:2562, 1993),
with rapid recovery from inactivation for potentials below rest.

Implementation by Paul B. Manis, April and Sept, 1999.
ENDCOMMENT

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(nA) = (nanoamp)
}

? interface
NEURON {
	SUFFIX jsrnaf
	USEION na READ ena WRITE ina
	RANGE gnac, gna, vsna
	RANGE minf, hinf, mtau, htau
	RANGE alpha_h, beta_h, alpha_m, beta_m
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	celsius = 22 (degC)
	tenC = 10 (degC)
	dt (ms)
	ena = 55 (mV)
	gnac = 0.07958 (mho/cm2)    <0,1e9>
	vsna = 0 (mV)
	q10 (1)
	qten(1)
	qtwo(1)}

STATE {
	m h
}

ASSIGNED {
	gna (mho/cm2)
    ina (mA/cm2)
    minf (1)
    hinf (1)
    mtau (ms) htau (ms)

	alpha_h(1)
	beta_h(1)
	alpha_m(1)
	beta_m(1)
}

LOCAL tablesave, mexp, hexp

INITIAL {
    rates(v)
    m = minf
    h = hinf

}

BREAKPOINT {
	SOLVE states
    gna = gnac*m*m*m*h
    ina = gna*(v - ena)
}

PROCEDURE states() {  :Computes state variables m, h, and n
	rates(v)      :             at the current v and dt.
	m = m + mexp*(minf-m)
	h = h + hexp*(hinf-h)
}


PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
LOCAL  x, tinc
	qtwo = 2^((celsius -22)/tenC)
	q10 = 3^((celsius - 22)/tenC)
	qten = 10^((celsius - 22)/tenC)

:"m" sodium activation system - JSR
	x = v + 49 + vsna
	alpha_m = 0.36 * q10 * vtrap(x, -x, 3.0(mV))
	x = v + 58 + vsna
	beta_m =  -0.4 * q10 * vtrap(x, x, 20(mV))
	mtau = 1(ms)/(alpha_m + beta_m)
	minf = alpha_m/(alpha_m + beta_m)

:"h" sodium inactivation system - JSR
	alpha_h = (2.4*q10/(1 + exp((v+68+vsna) / 3 (mV))))
	alpha_h = alpha_h + (0.8 * qten/( 1 + exp((v + 61.3 + vsna)/1 (mV))))
	beta_h = 3.6 * q10 / (1.0 + exp(-(v + 21 + vsna)/ 10 (mV)))
	htau = 1(ms)/(alpha_h + beta_h)
	hinf = alpha_h/(alpha_h + beta_h)

	tinc = -dt * q10
	mexp = 1 - exp(tinc/mtau)
	hexp = 1 - exp(tinc/htau)
}
UNITSOFF

FUNCTION vtrap(x (mV),y (mV), z(mV)) {  :Traps for 0 in denominator of rate eqns.
	if (fabs(y/z) < 1e-6) {
		vtrap = z*(1 - x/z/2)
	}else{
		vtrap = x/(1-exp(y/z))
	}
}
UNITSON

