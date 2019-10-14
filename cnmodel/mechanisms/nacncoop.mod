TITLE nacn.mod  A sodium conductance for a ventral cochlear nucleus neuron model

COMMENT

NEURON implementation of Jason Rothman's measurements of VCN conductances.

This file implements a modified version of the average brain sodium current
    used in the Rothman and Manis 2003 models.
    
The model differs from the one used in Rothman et al, (1993) in that the steep
voltage dependence of recovery from inactivation in that model is missing. This
may affect the refractory period. To use the other model, use jsrnaf.mod instead.

Original implementation by Paul B. Manis, April 1999 (JHU) and Sept 1999 (UNC-CH).

File split implementation, April 1, 2004.


Version nacncoop implements a cooperative sodium channel model built on the kinetics
of the original nacn model (R&M2003c). The motivation is to make a sodium channel with
faster activation kinetics, by introducing cooperativity between a subset of channels.
The model is based on concepts and implementation similar to Oz et al. 
J.Comp. Neurosci. 39: 63, 2015, and Huang et al., PloSOne 7:e37729, 2012.
The cooperative channels are modeled with the same kinetics as the non-cooperative
channels, but are treated as a separate subset (fraction: p). The cooperativity is
introduced by shifting the voltage "seen" by the channels by KJ*m^3*h, which moves
the channels to a faster regime (essentially, they experience a depolarized membrane
potential that depends on their current gating state, relative to the main population
of channels).

A subpopulation of Na channels (p [0..1]) experiences a small voltage-dependent shift
in the gating kinetics. The shift is determined by KJ

This version does not have all the temperature scaling. Does not pass modlunit.
Should work at 22C, appears to work at other temperatures ok.

Contact: pmanis@med.unc.edu

ENDCOMMENT

UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (nA) = (nanoamp)
}

NEURON {
        THREADSAFE
        SUFFIX nacncoop
        USEION na READ ena WRITE ina
        RANGE gbar, gna, ina, p, KJ
    	RANGE vshift : voltage shift parameter
        GLOBAL hinf, minf, htau, mtau, hinf2, minf2, htau2, mtau2
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
        v (mV)
        celsius (degC) : 22 (degC) model is defined at room temp in Baltimore
        dt (ms)
        ena (mV)
        gbar =  0.07958 (mho/cm2) <0,1e9>
        q10 = 3.0 : q10 for rates
        p  = 0.0 (): fraction of cooperative channels (0-1)
        KJ = 0 (mV) : coupling strength between cooperative channels (0-1000mV is usable range)
                      : setting either KJ = 0 or p = 0 will remove cooperativity.
        vshift = 0. (mV)
}

STATE {
        m h m2 h2
}

ASSIGNED {
    ina (mA/cm2) 
    gna (mho/cm2)
    vNa (mV)  : shifted V for cooperative behavior
    minf hinf minf2 hinf2
    mtau (ms) htau (ms) mtau2 (ms) htau2 (ms)

    }

LOCAL mexp, hexp, mexp2, hexp2

BREAKPOINT {
    SOLVE states METHOD cnexp
    
    gna = gbar*(p*(m2^3*h2) + (1.-p)*(m^3)*h)
    ina = gna*(v - ena)
}

UNITSOFF

INITIAL {
    rates(v)
    m = minf
    h = hinf
    m2 = minf2
    h2 = hinf2
    vNa = v + vshift + KJ*m^3*h
}

DERIVATIVE states {  :Computes state variables m, h, and n
    rates(v)      :             at the current v and dt.
    m' = (minf - m)/mtau 
    h' = (hinf - h)/htau
    m2' = (minf2 - m2)/mtau2
    h2' = (hinf2 - h2)/htau2
    vNa = v + vshift + KJ*m^3*h  : note addition of vshift here so that we do not add it in rates
}

LOCAL qt

PROCEDURE rates(v) {  :Computes rate and other constants at current v.
                      :Call once from HOC to initialize inf at resting v.

    qt = q10^((celsius - 22)/10) : if you don't like room temp, it can be changed!

: average sodium channel
    minf = 1 / (1+exp(-(v + 38 + vshift) / 7))
    hinf = 1 / (1+exp((v + 65 + vshift) / 6))
    mtau =  (10 / (5*exp((v + 60 + vshift) / 18) + 36*exp(-(v + 60+vshift) / 25))) + 0.04
    mtau = mtau/qt
    htau =  (100 / (7*exp((v + 60 + vshift) / 11) + 10*exp(-(v + 60 + vshift) / 25))) + 0.6
    htau = htau/qt

: cooperative group of channels
    minf2 = 1 / (1+exp(-(vNa + 38) / 7))
    hinf2 = 1 / (1+exp((vNa + 65) / 6))
    mtau2 =  (10 / (5*exp((vNa+60) / 18) + 36*exp(-(vNa+60) / 25))) + 0.04
    mtau2 = mtau2/qt
    htau2 =  (100 / (7*exp((vNa+60) / 11) + 10*exp(-(vNa+60) / 25))) + 0.6
    htau2 = htau2/qt

}

UNITSON
