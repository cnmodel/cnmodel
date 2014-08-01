TITLE KCNQ potassium channel for GPe neuron

COMMENT
 modeled by Gunay et al., 2008
 implemented in NEURON by Kitano, 2011
Threadsafe and unit checking, P.B. Manis, 2014

 ENDCOMMENT

UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
}

NEURON {
 THREADSAFE
 SUFFIX kcnq
 USEION k READ ek WRITE ik
 RANGE gbar, gk, iKCNQ
}

PARAMETER {
    v (mV)
    dt (ms)
    gbar  = 0.001 (mho/cm2)
    iKCNQ  = 0.0 (mA/cm2)
    ek (mV)

    theta_m = -61.0 (mV)
    k_m = 19.5 (mV)
    tau_m0 = 6.7 (ms)
    tau_m1 = 100.0 (ms)
    phi_m = -61.0 (mV)
    sigma_m0 = 35.0 (mV)
    sigma_m1 = -25.0 (mV)
}

STATE {
    m
}

ASSIGNED {
    ik (mA/cm2)
    gk (mho/cm2)
    minf
    taum (ms)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gk  = gbar*m*m*m*m
    ik = gk * (v-ek)
    iKCNQ = ik
}


INITIAL {
    settables(v)
    m = minf
}

DERIVATIVE states {
    settables(v)
    m' = (minf - m)/taum
}

PROCEDURE settables(v (mV)) {
    TABLE minf, taum FROM -100 TO 100 WITH 400

	minf = 1.0 / (1.0 + exp((theta_m - v)/k_m))
	taum = tau_m0 + (tau_m1 - tau_m0)/(exp((phi_m - v)/sigma_m0) + exp((phi_m - v)/sigma_m1))
}
