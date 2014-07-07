TITLE Low threshold calcium current
:
:   Ca++ current responsible for low threshold spikes (LTS)
:   RETICULAR THALAMUS
:   Differential equations
:
:   Model of Huguenard & McCormick, J Neurophysiol 68: 1373-1383, 1992.
:   The kinetics is described by standard equations (NOT GHK)
:   using a m2h format, according to the voltage-clamp data
:   (whole cell patch clamp) of Huguenard & Prince, J Neurosci.
:   12: 3804-3817, 1992.
:
:    - Kinetics adapted to fit the T-channel of reticular neuron
:    - Time constant tau_h refitted from experimental data
:    - shift parameter for screening charge
:
:   Model described in detail in:   
:     Destexhe, A., Contreras, D., Steriade, M., Sejnowski, T.J. and
:     Huguenard, J.R.  In vivo, in vitro and computational analysis of
:     dendritic calcium currents in thalamic reticular neurons.
:     Journal of Neuroscience 16: 169-185, 1996.
:
:
:   Written by Alain Destexhe, Salk Institute, Sept 18, 1992
:

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX it2
	USEION ca READ cai, cao WRITE ica
	RANGE gcabar, m_inf, tau_m, h_inf, tau_h, shift, qm, qh
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)

	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v		(mV)
	celsius	= 36	(degC)
:	eca	= 120	(mV)
	gcabar	= .003	(mho/cm2)
	shift	= 0 	(mV)
	cai	= 2.4e-4 (mM)		: adjusted for eca=120 mV
	cao	= 2	(mM)
	qm	= 2.5
	qh 	= 2.5
}

STATE {
	m h
}

ASSIGNED {
	ica	(mA/cm2)
	carev	(mV)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	phi_m
	phi_h
}

BREAKPOINT {
	SOLVE castate METHOD euler
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcabar * m*m*h * (v-carev)
}

DERIVATIVE castate {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {
:
:   Activation functions and kinetics were obtained from
:   Huguenard & Prince, and were at 23-25 deg.
:   Transformation to 36 deg using Q10
:
	phi_m = qm ^ ((celsius-24)/10)
	phi_h = qh ^ ((celsius-24)/10)

	evaluate_fct(v)
	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { 
:
:   Time constants were obtained from J. Huguenard
:

	m_inf = 1.0 / ( 1 + exp(-(v+shift+50)/7.4) )
	h_inf = 1.0 / ( 1 + exp((v+shift+78)/5.0) )

	tau_m = ( 3 + 1.0 / ( exp((v+shift+25)/10) + exp(-(v+shift+100)/15) ) ) / phi_m
	tau_h = ( 85 + 1.0 / ( exp((v+shift+46)/4) + exp(-(v+shift+405)/50) ) ) / phi_h
}
UNITSON
