TITLE Moore-Cox sodium channel
: Biophy. J. (1976) 16:171
: This paper mapped HH VClamp currents, no attention paid to resting currents. 
: Used V-jump @ T=0 to check match with HH action  potential; fig shows equal V-threshold levels
: Ramon first noted instability at rest, spontaneous impulse generation.
: Same problem noted with resimulation with NEURON. Now thresholds rather different
: Revised July 28, 1995 to remove instability. Added back reaction rate coefficients for HH beta m
: First use of  NEURON's new "Run Fitter" to find best values of these coefficients,
: using delay and fitting  ina to just beyond  peak.
: Excellent fit to HH AP with these coefficients except for "gratituitous hump" in HH
: Changed from HH reference potential level at rest to NEURON @-65mV
: Added parameter vrest to allow adjustment of relative voltages for other simulations P. Manis 1/22/99

?interface
NEURON {
     SUFFIX MCna
     USEION na READ ena WRITE ina
     RANGE gna, ina, lp, ml, nm, porate     
	 GLOBAL gnabar
}

UNITS {
     (mA) = (milliamp)
     (mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
     v (mV)
     celsius = 22  (degC)
     dt (ms)
     gnabar=.120 (mho/cm2)
     ena = 50 (mV)
     lp=1.9
     ml=0.75
     nm=0.3
     porate=1
	 vrest = -55 (mV) : define subtraction point for rest
}
STATE {
     P L M N O
}
ASSIGNED {
     ina (mA/cm2)
     gna (mho/cm2)
}

LOCAL am, bm, ah, bh

INITIAL {
     P=1
     rate(v*1(/mV))
     SOLVE states STEADYSTATE sparse
}

BREAKPOINT {
     SOLVE states METHOD sparse
     ina = gnabar*N*(v - ena)
}

KINETIC states {
     rate(v*1(/mV))
     CONSERVE P + L + M + N + O = 1
     ~ P <-> L (am, lp*bm)    :back reaction in original = 3.5   
     ~ L <-> M (2*am, ml*bm)          :back reaction in original = 0
     ~ M <-> N (3*am, nm*bm)         :back reaction in original = 0
     ~ N <-> O (1.1*bh, 0)
     ~ N <-> P (3*bm, 0)
     ~ P <-> O (bh, ah)       :back reaction in original = 1,
  			      : found this to still be good
}

UNITSOFF
FUNCTION alp(v(mV),i) { LOCAL a,b,c,q10 : order m,h
     v = -v + vrest :convert to hh convention
     q10 = 3^((celsius - 6.3)/10)
     if (i==0) {
          alp = q10*0.1*expM1(v + 25, 10)
     }else if (i==1){
          alp = q10*0.07*exp(v/20)
     }
}

FUNCTION bet(v,i) { LOCAL a,b,c,q10 : order m,h
     v = -v + vrest
     q10 = 3^((celsius - 6.3)/10)
     if (i==0) {
          bet = q10* 4*exp(v/18)
     }else if (i==1){
          bet = q10*1/(exp(0.1*v + 3) + 1)
     }
}

FUNCTION expM1(x,y) {
     if (fabs(x/y) < 1e-6) {
          expM1 = y*(1 - x/y/2) : for singular point
     }else{
          expM1 = x/(exp(x/y) - 1)
     }
}

PROCEDURE rate(v) {LOCAL a, b, tau :
     TABLE am, ah, bm, bh DEPEND dt, celsius FROM -100 TO 100 WITH 200
     am = alp(v, 0)
     ah = alp(v, 1)
     bm = bet(v, 0)
     bh = bet(v, 1)
}
UNITSON






