: AdEx GIF model

: This implementation is for the equations in:
: Naud R, Marcille N, Clopath C, Gerstner W. Firing patterns in the adaptive
: exponential integrate-and-fire model. Biol Cybern. 2008 Nov;99(4-5):335-47. doi: 
: 10.1007/s00422-008-0264-7. Epub 2008 Nov 15. PubMed PMID: 19011922; PubMed
: Central PMCID: PMC2798047.

: Which in turn is based on:
: Brette R, Gerstner W. Adaptive exponential integrate-and-fire model as an
: effective description of neuronal activity. J Neurophysiol. 2005
: Nov;94(5):3637-42. Epub 2005 Jul 13. PubMed PMID: 16014787.
:
: Paul B. Manis
: 9 Nov 2017, Washington DC
:


NEURON {
:	ARTIFICIAL_CELL AdEx
    SUFFIX AdEx
    RANGE gl, el, delt, vt, vr, w, b, cm, is, a, tauw
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
    vt = -50 (mV)  
    vr = -58 (mV): reset value after a spike
    a = 2 (nS)
    b = 0 (pA)
    is = 0 (pA)
    tauw = 30 (ms)
	refract = 1 (ms)
}

ASSIGNED {
    i       (mA/cm2)
	t0(ms)  : time of last spike
	refractory  : flag indicating when in a refractory period
}

STATE {
    w
    Vm
}

INITIAL {
	Vm = el
	t0 = t
    w = 0
    
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
    w' = (a*(Vm - el) - w)/tauw
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
    Vm' = ( -gl*(Vm-el) + gl*delt*et + is - w)/cm
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