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
THREADSAFE
	SUFFIX jsrna
	USEION na READ ena WRITE ina
	RANGE gbar
    RANGE gna, vsna
	RANGE minf, hinf, mtau, htau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	v (mV)
	dt (ms)
    ena = 55 (mV)
	gbar = 0.25 (mho/cm2)    <0,1e9>
	vsna = 0 (mV)
    q10 = 3.0 (1)
}

STATE {
	m h 
}

ASSIGNED {
	gna (mho/cm2) 
    ina (mA/cm2)
    minf hinf 
    mtau (ms) htau (ms)
    celsius (degC)
}

LOCAL mexp, hexp 

? currents
BREAKPOINT {
	SOLVE states METHOD cnexp
    gna = gbar*(m^3)*h
    ina = gna*(v - ena)
}

UNITSOFF

INITIAL {
    trates(v)
    m = minf
    h = hinf
}

DERIVATIVE states {  :Computes state variables m, h, and n
	trates(v)      :             at the current v and dt.
	m' = (minf - m)/mtau
    h' = (hinf - h)/htau
}

LOCAL qt

PROCEDURE rates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
LOCAL  alpha, beta, sum

	qt = q10^((celsius - 22)/10 (degC)) : R&M'03 used 3

: Note qt temperature here cancels in minf (a/(a+b))
:"m" sodium activation system - JSR
        alpha = -0.36*qt*vtrap((v+49),-3)
        beta =  0.4*qt*vtrap((v+58),20)
        sum = alpha + beta
		mtau = 1/sum
        minf = alpha/sum

:"h" sodium inactivation system - JSR
        alpha = 2.4*qt/(1+exp((v+68-vsna)/3 (mV))) + 0.8*qt/(1+exp(v+61.3-vsna))
        beta = 3.6*qt/(1+exp(-(v+21-vsna)/10 (mV)))
        sum = alpha + beta
		htau = 1/sum
        hinf = alpha/sum




: jsr modified sodium channel - defined in terms of alpha and beta this time
:    am = (0.36*q10*(v+49))/(1-exp(-((v+49)/3)))
:    am = -(0.36*q10*vtrap(-(v+49),3))
:    bm = -(0.40*q10*(v+58))/(1-exp((v+58)/20))
:    bm = (0.40*q10*vtrap((v+58),20))
:    ah = ((2.4*q10)/(1+exp((v+68)/3))) + (0.8*qten/(1+exp(v+61.3)))
:    bh = (3.6*q10)/(1+exp(-(v+21)/10))

:    minf = am/(am+bm)
:    hinf = ah/(ah+bh)
    
:	mtau = 1/(am+bm)
:	htau = 1/(ah+bh)

}

PROCEDURE trates(v (mV)) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.
	LOCAL tinc

    rates(v)    : not consistently executed from here if usetable_hh == 1
        : so don't expect the tau values to be tracking along with
        : the inf values in hoc

	tinc = -dt : * q10 # handled in rates now
	mexp = 1 - exp(tinc/mtau)
	hexp = 1 - exp(tinc/htau)
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(exp(x/y) - 1)
	}
}


