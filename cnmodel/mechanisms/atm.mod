: ATM GIF model

: This implementation the adaptive theshold model (ATM) is for the equations in:
: Fontaine, B., Benichourx, V., Joris, P.X., and Brette, R. Prediciting
: spike timing in hhigy synchronous auditory neurons at different sound
: levels. J. Neurophysiol. 110: 1672-1688, 2013.

: Which in turn is based on:
: Brette R, Gerstner W. Adaptive exponential integrate-and-fire model as an
: effective description of neuronal activity. J Neurophysiol. 2005
: Nov;94(5):3637-42. Epub 2005 Jul 13. PubMed PMID: 16014787.
:
: Paul B. Manis
: 2 December 2017, Chapel Hill, NC
:
: Incomplete version

NEURON {
:	ARTIFICIAL_CELL ATM
    SUFFIX ATM
    RANGE gl, el, delt, vt, vr, alpha, beta, cm, is, a, tauw
    RANGE refract, Vm
    NONSPECIFIC_CURRENT i
	: m plays the role of voltage
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	cm = 200 (pF)
    el = -70 (mV)  : leak (RMP)
    gl = 10 (nS)  : resting input R
    delt = 2 (mV)  : spike threshold sharpness
    vr = -58 (mV): reset value after a spike
    a = 2 (1)
    b = 2 (1)
    beta = 0 (1)
    alpha = 0
    is = 0 (pA)
    taut = 30 (ms) : threshold tau
	refract = 1 (ms)
}

ASSIGNED {
    i       (mA/cm2)
	t0      (ms)    : time of last spike
	refractory      : flag indicating when in a refractory period
}

STATE {
    w
    Vm
    vt
}

INITIAL {
	Vm = el
	t0 = t
    a = 0
    b = 0
    
	refractory = 0 : 0-integrates input, 1-refractory
}

BREAKPOINT {
    SOLVE states METHOD cnexp

    if (refractory == 0 && Vm <= 0.) {
        states()
    }
    if (refractory == 1) {
        if ((t-t0) >= refract){
            refractory = 0
            Vm = vr
            states()
        }
        else {
            Vm = 0.
        }
    }
    if (refractory == 0 && Vm > 0.) {
        refractory = 1
        t0 = t
        Vm = 0.
        w = w + b
    }
}


DERIVATIVE states {  : update adaptation variable w
    LOCAL eterm, et
    vt' = (a*i - vt)/taut
COMMENT
    eterm = (Vm-vt)/delt
    if (eterm > 700 ) { : prevent overflow of the exponential term 
                        : (it would be better to estimate the value... but for this 
                        : implementation, not necessary as this will be the term
                        : that drives the model to spike - after that V is reset
                        : so the time evolution no longer matters)
         et = 700.
    }
    else {
        et = exp(eterm)
    }
ENDCOMMENT
    Vm' = gl*( -(Vm-el) + i)/cm
}


COMMENT
NET_RECEIVE (w) {
	if (refractory == 0) { : inputs integrated only when excitable
        i = -gl*(v-el) + gl*delt*exp((Vm-vt)/delt) - w
        m = i/cm
		t0 = t
        states()
        if (m > 0) {
			refractory = 1
			m = 0
			net_send(refractory, refractory)
			net_event(t)
		}
	} else if (flag == 1) { : ready to integrate again
		t0 = t
		refractory = 0
		m = vr
	}
}
ENDCOMMENT