TITLE ihsgc.mod - Spiral Ganglion Cell Ih current

COMMENT
Ih for Spiral ganglion cells.
Kinetcs are based on average fits to mouse SGCs, all regions
(apical, middle and basal).
Data used to establish the kinectic parameters were collected by
Qing Liu and Robin Davis (Rutgers).
Data were taken at room temperature.
Kinetic parameters were extracted by curve fitting for fast and
slow components from activation and deactivation (using
the program Ihfit4.py).

Implementation by Paul B. Manis, January-April, 2012.
pmanis@med.unc.edu

Note: vshift parameter is nominally 0. This parameter can
shift the entire activation and rate curves, keeping them
in register for each component of the conductance.

ENDCOMMENT

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(nA) = (nanoamp)
}

NEURON {
	SUFFIX ihsgc
	NONSPECIFIC_CURRENT i
	RANGE ghbar, gh, ih, eh, vshift
	RANGE vh, k, vhs, ks
	RANGE rinf, rtau, sinf, stau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
		v (mV)
		celsius = 22 (degC)
		dt (ms)
		ghbar = 0.00318 (mho/cm2) <0,1e9>

: Parameters from kinetic analysis
: Format for NEURON MOD file: 

         : (Run on date =  2013-10-17 17:28:10.659747 )

: A. Fast component (Fast trace):
        : These parameters could be automatically replaced with the new fit values
        : from the analysis. 
: D =   29.33, A1= 0.017486, v1 =   61.25, k1 = 26.2780, A2 = 0.000027, v2 = 30.7333, k2 = 11.2766
                : Boltzmann:
                eh = -43 (mV)
                vh = -103.035 (mV)
                k =  11.223 (mV)
                vshift = 0.0 (mV)
                afast =  0.5350 : fractiom that is fast.

                :tau
                taufac = 1.0 (1)
                taumin = 29.3334 (ms)
                tausc1 = 0.0174863 (ms)
                vtau1 = 61.2494 (mV)
                kfac1 = 26.278 (mV)
                tausc2 = 2.65845e-05 (ms)
                vtau2 = 30.7333 (mV)
                kfac2 = 11.2766 (mV)

        
: B. Slow component (Cyan trace):
            : (Run on date =  2013-10-17 17:28:10.659849 )

: D =   40.78, A1= 0.007921, v1 =   30.74, k1 = 37.3538, A2 = 0.000011, v2 = 13.9087, k2 = 19.3753
                : Boltzmann:
                : eh = -43 (mV)
                svh1 = -87.981 (mV)
                sk1 =   2.343 (mV) : double boltzmann 
                svh2 = -116.164 (mV)
                sk2 =  -5.281 (mV)
                svshift = 0.0 (mV)

                sba1 = 0.704464 : amplitdue first component (0 .. 1)
                sba2 =0.295536  : amplitude slow component (1 - sba2)
                aslow =  0.4650 : total slow 
                :stau
                staufac = 1.0 (1)
                staumin = 40.7778 (ms)
                stausc1 = 0.00792101 (ms)
                svtau1 = 30.7426 (mV)
                skfac1 = 37.3538 (mV)
                stausc2 = 1.12963e-05 (ms)
                svtau2 = 13.9087 (mV)
                skfac2 = 19.3753 (mV)

}

STATE {
		r
		s
}

ASSIGNED {
	gh (mho/cm2)
	i (mA/cm2)
	ih (mA/cm2)
	rinf
	rtau (ms)
	sinf
	stau (ms)
}

LOCAL rexp, sexp

BREAKPOINT {
	SOLVE states
	gh = ghbar*(afast*(r^2)+aslow*s) : Balance between fast and slow determined by afast and aslow
	ih = gh*(v - eh)
	i = ih
}

UNITSOFF

INITIAL {
	trates(v)
	r = rinf
	s = sinf
}

PROCEDURE states() {  : Updates state variables r and s
	trates(v)		  : at the current voltage
	r = r + rexp*(rinf-r)
	s = s + sexp*(sinf-s)
VERBATIM
	return 0;
ENDVERBATIM
}

LOCAL q10, st
PROCEDURE rates(v) {  : Computes rate and activation at voltage = v.

	q10 = 3.0^((celsius - 22.0)/10.0) : adjust for temperature...

: fast component - standard HH-like kinetics.
	rinf = 1.0 / (1+exp((v - vh + vshift) / k))^0.5
	rtau = tausc1*exp((v + vtau1 + vshift) / kfac1) + tausc2*exp(-(v + vtau2 + vshift) / kfac2)
	rtau = (taumin + taufac/rtau)

: slow component
: double boltzman activation function (decreasing conductance), unequal sharing. 
	sinf = 1. / (1+exp((v - svh1 + vshift) / sk1))
	st = 1. / (1 + exp((v - svh2 + vshift) / sk2))
	sinf = sinf + sba2*(st - 1) : normalize
	stau = staufac / (stausc1*exp((v + svtau1 + vshift) / skfac1) + stausc2*exp(-(v + svtau2 + vshift) / skfac2))
	stau = (stau + staumin)
}

PROCEDURE trates(v) {	: Computes rate and other constants at voltage v.
	LOCAL tinc
	TABLE rinf, rexp, sinf, sexp
	DEPEND dt, celsius FROM -200 TO 150 WITH 350

	rates(v)
	tinc = -dt * q10
	rexp = 1.0 - exp(tinc/rtau)
	sexp = 1.0 - exp(tinc/stau)
}

UNITSON
