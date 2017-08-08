TITLE h current for Octopus cells of Cochlear Nucleus
: From Bal and Oertel (2000)

: Modified, P. Manis July 2014, 2017 
: Parameters from McGinley et al. paper

NEURON {
    THREADSAFE
    SUFFIX hcnobo
    NONSPECIFIC_CURRENT i
    RANGE  gbar, eh, gh, q10tau
    GLOBAL hinf, tau1, tau2
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (um) = (micron)
    R = (k-mole)(joule/degC)
    F = (faraday)(kilocoulombs)
}

PARAMETER {
    gbar = 0.0005       (mho/cm2)   
                                
    vhalf1  = -70 (mV) :  -50   (mV)        : v 1/2 for forward
    vhalf2  = -84   (mV)        : v 1/2 for backward    
    gm1   = 0.3 :(mV)           : slope for forward
    gm2   = 0.6  :    (mV)      : slope for backward
    zeta1   = 3 :   (/ms)       
:    zeta2   = 3 :   (/ms)       
    a01 = 4.8e-3 (/ms)  : was 0.008
    a02 = 2.9e-3 (/ms) : was 0.0029 (/ms)
    frac = 0.8
    c0 = 273.16  (degC)
    thinf  = -72.4    (mV)        : inact inf slope   
    qinf  = 5.3   (mV)        : inact inf slope 
    q10tau = 4.5             : from Magee (1998)
    v       (mV)
}


ASSIGNED {
    celsius (degC)
    i       (mA/cm2)
    gh      (mho/cm2)
    eh      (mV)    : must be explicitly def. in hoc
    hinf 
    tau1 (ms)
    tau2 (ms) 
    q10 ()
	ssih
    ct
}
 

STATE { h1 h2 }

BREAKPOINT {
    SOLVE states METHOD cnexp
    : SOLVE states METHOD derivimplicit
    gh = gbar*(h1*frac + h2*(1.0-frac))
    i = gh * (v - eh)
} 

INITIAL {
    ct = 1e-3*zeta1*F/(R*(c0+celsius))
    
    q10 = q10tau^((celsius - 33.0)/10.0 (degC)) : Measurements at 33
    rates(v)
    h1=hinf
    h2=hinf
	ssih = 0.
}

DERIVATIVE states {   
    rates(v)      
    h1' = (hinf - h1)/tau1
    h2' = (hinf - h2)/tau2
}

PROCEDURE rates(v (mV)) {  
    tau1 = bet1(v)/(q10*a01*(1.0+alp1(v)))
    tau2 = bet2(v)/(q10*a02*(1.0+alp2(v)))
    hinf = 1.0/(1.0+exp((v-thinf)/qinf))
}

FUNCTION alp1(v(mV)) {
  alp1 = exp((v-vhalf1)*ct) 
}

FUNCTION bet1(v(mV)) {
  bet1 = exp(gm1*(v-vhalf1)*ct)
}

FUNCTION alp2(v(mV)) {
  alp2 = exp((v-vhalf2)*ct) 
}

FUNCTION bet2(v(mV)) {
  bet2 = exp(gm2*(v-vhalf2)*ct) 
}
