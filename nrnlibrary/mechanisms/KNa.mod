COMMENT
-----------------------------------------------------------------------------
From "Slack and Slick KNa Channels Regulate the Accuracy of Timing of Auditory Neurons"
Bo Yang, Rooma Desai, and Leonard K. Kaczmarek
J Neuroscience 27(10)2617-2627 (2007)
and
Localization of a High Threshold Potassium Channel in the Rat Cochlear Nucleus
TERESA M. PERNEY1* AND LEONARD K. KACZMAREK2
J Comp Neurol 386:178-202 (1997)


-----------------------------------------------------------------------------
ENDCOMMENT


TITLE kna.mod

NEURON {
	SUFFIX kna
	USEION k READ ek WRITE ik
	USEION na READ nai
	RANGE gkna, ninf, alpha, beta, sinf, stau
	RANGE etaalpha, etabeta, kalpha, kbeta, kfs, kbs
	GLOBAL gknabar
}


PARAMETER {
	gknabar = .0000010	(mho/cm2)
	:gknabar = 0	(mho/cm2)  
	v 		(mV)
	nai = 2	(mM)
	etaalpha = 0.0558 (/mV)
	etabeta = 0.0308 (/mV)
	kalpha = 0.141 (/ms)
	kbeta = 0.0587 (/ms)
	kfs = 0.14 (/mM ms)
	kbs = 2.0 (/mM ms)
	ek = -91 (mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) =	(milli/liter)
} 

ASSIGNED {
	ina 		(mA/cm2)
	gkna		(mho/cm2)
	ninf
	sinf
	alpha
	beta
	stau
	ik (mA/cm2)
}
 

STATE {n s}

UNITSOFF

BREAKPOINT {
    SOLVE states METHOD cnexp
	gkna =gknabar*n*n*s*s*s*s
	ik = gkna*(v+80)
}
	
DERIVATIVE states {
	rates(v,nai)
	n' = (alpha*(1-ninf))-beta
	s' = (sinf - s)/stau
}	
	

INITIAL {
	rates(v,nai)
	s = sinf
	n = ninf	
}

PROCEDURE rates(v(mV),nai(mM)) {
	alpha = knalpha(v)
	beta = knabeta(v)
	ninf = Kna_n(alpha,beta)
	sinf = Kna_s(nai)
	stau = tau(nai)
}

FUNCTION Kna_s(nai (mM)) {  :Sodium dependence variable
	Kna_s = (kfs*nai)/(kbs + (kfs*nai))
}

FUNCTION tau(nai (mM)){
	tau = 1/(kbs + (kfs*nai))
}

FUNCTION Kna_n(alpha,beta) {  :Gating variable
	Kna_n = alpha/(alpha + beta)
}

FUNCTION knalpha(x) {  :Gating variable
	knalpha = kalpha*(exp(etaalpha*x))
}

FUNCTION knabeta(x) {  :Gating variable
	knabeta = kbeta*(exp(etabeta*x))
}


UNITSON
