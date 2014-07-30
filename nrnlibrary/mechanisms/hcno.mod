TITLE h current for Octopus cells of Cochlear Nucleus
: From Bal and Oertel (2000)
: M.Migliore Oct. 2001
: Modified, P. Manis July 2014.

NEURON {
    THREADSAFE
    SUFFIX hcno
    NONSPECIFIC_CURRENT i
    RANGE  gbar, eh, gh
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
                                
    vhalf1  = -50   (mV)        : v 1/2 for forward
    vhalf2  = -84   (mV)        : v 1/2 for backward    
    gm1   = 0.3 :(mV)           : slope for forward
    gm2   = 0.6  :    (mV)      : slope for backward
    zeta1   = 3 :   (/ms)       
    zeta2   = 3 :   (/ms)       
    a01 = 0.008  (/ms)
    a02 = 0.0029 (/ms)
    frac=0.0
    c0 = 273.16  (degC)
    thinf  = -66    (mV)        : inact inf slope   
    qinf  = 7   (mV)        : inact inf slope 
    q10tau = 4.5             : from Magee (1998)
    v       (mV)
    q10g = 2.0    : Rothman...
}


ASSIGNED {
    celsius (degC)
    i       (mA/cm2)
    gh      (mho/cm2)
    eh      (mV)    : must be explicitly def. in hoc
    hinf 
    tau1 (ms)
    tau2 (ms) 
    qg ()  : computed q10 for gnabar based on q10g
    q10 ()
	ssih
}
 

STATE { h1 h2 }

BREAKPOINT {
    SOLVE states METHOD derivimplicit
    gh = qg*gbar*(h1*frac + h2*(1-frac))
    i = gh * (v - eh)
} 

INITIAL {
    qg = q10g^((celsius-33)/10 (degC))  :note original measurements made at 33 C
    q10 = q10tau^((celsius - 22)/10 (degC)) : if you don't like room temp, it can be changed!
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
    tau1 = bet1(v)/(q10*a01*(1+alp1(v)))
    tau2 = bet2(v)/(q10*a02*(1+alp2(v)))
    hinf = 1/(1+exp((v-thinf)/qinf))
}

FUNCTION alp1(v(mV)) {
  alp1 = exp(1e-3*zeta1*(v-vhalf1)*F/(R*(c0+celsius))) 
}

FUNCTION bet1(v(mV)) {
  bet1 = exp(1.e-3*zeta1*gm1*(v-vhalf1)*F/(R*(c0+celsius)))
}

FUNCTION alp2(v(mV)) {
  alp2 = exp(1.e-3*zeta2*(v-vhalf2)*F/(R*(c0+celsius))) 
}

FUNCTION bet2(v(mV)) {
  bet2 = exp(1.e-3*zeta2*gm2*(v-vhalf2)*F/(R*(c0+celsius))) 
}
