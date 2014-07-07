COMMENT
This T-type calcium current was originally reported in Wang XJ et al 1991
This file supplies a version of this current identical to Quadroni and Knopfel 1994
except for gbar and Erev (see notes below).
ENDCOMMENT

NEURON {
	SUFFIX lva
	: NONSPECIFIC_CURRENT i
	USEION ca WRITE  ica
	RANGE Erev,g, gbar, i
	RANGE k, taum, minf, alpha_1, alpha_2, beta_1, beta_2, V_s
}

UNITS {
	(S)	=	(siemens)
	(mV)	=	(millivolt)
	(mA)	=	(milliamp)
}

PARAMETER {
	gbar = 0.4e-3	(S/cm2) < 0, 1e9 > : Quadroni and Knopfel use 166e-6
	Erev = 120 (mV)	: orig from Wang XJ et al 1991 was 120
	: note: Quadroni and Knopfel 1994 table 1 use 80 instead
	V_s = 0 (mV)	: used to describe effect of changing extracellular [Ca]
			: 0 corresponds to [Ca]outside = 3 mM (p 841)
}

ASSIGNED {
	ica (mA/cm2)
	i (mA/cm2)
	v (mV)
	g (S/cm2)
	k
	alpha_1 (1)
	alpha_2	(1)
	beta_1 (1)
	beta_2 (1)
}

STATE {	m h d }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar * m^3 * h
	ica = g * (v - Erev)
	i = ica	: used only to display the value of the current (section.i_lva(0.5))
}

INITIAL {
	LOCAL C, E
	: assume that v has been constant for a long time
	: (derivable from rate equations in DERIVATIVE block at equilibrium)
	rates(v)
	m = minf(v)
	: h and d are intertwined so more complex than above equilib state for m
	C =  beta_1 / alpha_1
	E =  alpha_2 / beta_2
	h = E / (E * C + E + C)
	d = 1 - (1 + C) * h
}

DERIVATIVE states{ 
	rates(v)
	m' = (minf(v) - m)/taum(v)		: alpham(v) * (1 - m) - betam(v) * m
	h' = alpha_1 * (1 - h - d) - beta_1 * h
	d' =  beta_2 * (1 - h - d) - alpha_2 * d
}

FUNCTION minf(Vm (mV)) (1) {
	minf = 1.0 / (1.0 + exp(-(Vm + V_s + 63)/7.8))
}

FUNCTION taum(Vm (mV)) (ms) {
	taum = (1.7 + exp( -(Vm + V_s + 28.8)/13.5 )) / (1.0 + exp( -(Vm + V_s + 63)/7.8) )
}

PROCEDURE rates(Vm(mV)) { LOCAL tau_2
	k = (0.25 + exp((Vm + V_s + 83.5)/6.3))^0.5 - 0.5
	tau_2 = 240.0 / (1 + exp((Vm + V_s + 37.4)/30)) : same as tau2 p 842 equation (15)
	alpha_1 = exp( -(Vm + V_s +160.3)/17.8 )	: p 842  equation (14)
	beta_1 = k * alpha_1
	alpha_2 = 1.0 / ( tau_2 * (1.0 + k) )
	beta_2 = k * alpha_2
}