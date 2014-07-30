TITLE Slow Ca-dependent potassium current
:
:   Ca++ dependent K+ current IC responsible for slow AHP
:   Differential equations
:
:   Model based on a first order kinetic scheme
:
:      <closed> + n cai <-> <open>	(alpha,beta)
:
:   Following this model, the activation fct will be half-activated at 
:   a concentration of Cai = (beta/alpha)^(1/n) = cac (parameter)
:
:   The mod file is here written for the case n=2 (2 binding sites)
:   ---------------------------------------------
:
:   This current models the "slow" IK[Ca] (IAHP): 
:      - potassium current
:      - activated by intracellular calcium
:      - NOT voltage dependent
:
:   A minimal value for the time constant has been added
:
:   Ref: Destexhe et al., J. Neurophysiology 72: 803-818, 1994.
:
:   Modifications by Arthur Houweling for use in MyFirstNEURON


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	THREADSAFE
    SUFFIX kpksk
	USEION k READ ek WRITE ik
	USEION ca READ cai
    RANGE  m_inf, tau_m, gbar, gk
	GLOBAL beta, cac
	RANGE ik
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}


PARAMETER {
	v		(mV)
	celsius		(degC)
    dt      (ms)
	ek		(mV)
	cai		(mM)	
	gbar	= .01	(mho/cm2)
	beta	= 0.002	(1/ms)		: backward rate constant
	cac	= 0.010	(mM)		: middle point of activation fct
	taumin	= 0.1	(ms)		: minimal value of the time cst
}


STATE {
	m
}

ASSIGNED {
	ik	(mA/cm2)
    gk  (mho/cm2)
	m_inf
	tau_m	(ms)
	tadj ()
}

BREAKPOINT { 
	SOLVE states :METHOD euler
	gk = gbar * m*m
    ik =  gk * (v - ek)
}

:DERIVATIVE states {
:       evaluate_fct(v,cai)
:
:       m'= (m_inf-m) / tau_m 
:}
  
PROCEDURE states() {
        evaluate_fct(v,cai)

        m= m + (1-exp(-dt/tau_m))*(m_inf-m)
}

INITIAL {
:
:  activation kinetics are assumed to be at 22 deg. C
:  Q10 is assumed to be 3
:
	tadj = 3 ^ ((celsius-22.0)/10 (degC))

	evaluate_fct(v,cai)
	m = m_inf
}

PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL car

	car = (cai/cac)^2

	m_inf = car / ( 1 + car )
	tau_m = 1 / beta / (1 + car) / tadj

    if(tau_m < taumin) { tau_m = taumin } 	: min value of time cst
}
