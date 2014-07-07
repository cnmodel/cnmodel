TITLE hyperpolarization-activated current (H-current)

COMMENT
	"Standard iH: ih_std"
	Two distinct activation gates are assumed with the same asymptotic
	opening values, a fast gate (F) and a slow gate (S). The following
	kinetic scheme is assumed

	s0  --(Alpha)--> s1 + n Cai  --(k1)--> s2
           <--(Beta)---             <--(k2)--

 	f0  --(Alpha)--> f1 + n Cai  --(k1)--> f2
           <--(Beta)---             <--(k2)--

	where s0/f0, s1/f1, and s2/f2 are resp. fraction of closed slow/fast
	gates, fraction of open unbound slow/fast gates, and fraction of open
	calcium-bound slow/fast	gates, n is taken 2, and k1 = k2*C where
	C = (cai/cac)^n and cac is the critical value at which Ca2+ binding
	is half-activated.

	The total current is computed according

	ih = ghbar * (s1+s2) * (f1+f2) * (v-eh)

        *********************************************
        reference:      Destexhe, Babloyantz & Sejnowski (1993)
			Biophys.J. 65, 1538-1552
        found in:       thalamocortical neurons
        *********************************************
	Maxim Bazhenov's first mod file
        Rewritten for MyFirstNEURON by Arthur Houweling
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iH_std
	USEION h READ eh WRITE ih VALENCE 1
	USEION ca READ cai
    RANGE ghbar, tau_s, tau_f, tau_c, ih_std, eh, vshift
	GLOBAL cac
}

UNITS {
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
	(molar)	= (1/liter)
	(mM) 	= (millimolar)
}

PARAMETER {
	v		    (mV)
	cai		    (mM)
	celsius		(degC)
	eh	= -43	(mV)
	ghbar = 4e-5	(mho/cm2)
	cac  = 5e-4	(mM)
	vshift = 0 (mV)
}

STATE {
	s1			: fraction of open unbound slow gates
	s2 			: fraction of open calcium-bound slow gates
	f1	 		: fraction of open unbound fast gates
	f2			: fraction of open calcium-bound fast gates
}

ASSIGNED {
    ih		    (mA/cm2)
    h_inf
    tau_s		(ms)	: time constant slow gate
    tau_f 		(ms)	: time constant fast gate
    tau_c		(ms)	: time constant calcium binding
    alpha_s		(1/ms)
    alpha_f 	(1/ms)
    beta_s 		(1/ms)
    beta_f		(1/ms)
    C
    k2		    (1/ms)
    tadj
    s0			: fraction of closed slow gates
    f0			: fraction of closed fast gates
}

BREAKPOINT {
	SOLVE states METHOD euler
	ih = ghbar * (s1+s2) * (f1+f2) * (v-eh)
}

UNITSOFF
DERIVATIVE states {
	evaluate_fct(v,cai)

    s1' = alpha_s*s0 - beta_s*s1 + k2*(s2-C*s1)
    f1' = alpha_f*f0 - beta_f*f1 + k2*(f2-C*f1)
    s2' = -k2*(s2-C*s1)
    f2' = -k2*(f2-C*f1)

    s0 = 1-s1-s2
    f0 = 1-f1-f2
}

INITIAL {
	: Q10 assumed to be 3
	tadj = 3^((celsius-35.5)/10)
	evaluate_fct(v,cai)

	s1 = alpha_s / (beta_s+alpha_s*(1+C))
	s2 = alpha_s*C / (beta_s+alpha_s*(1+C))
	s0 = 1-s1-s2
	f1 = alpha_f / (beta_f+alpha_f*(1+C))
	f2 = alpha_f*C / (beta_f+alpha_f*(1+C))
	f0 = 1-f1-f2

	tau_c = 1 / (1+C) / k2	: for plotting purposes
}

PROCEDURE evaluate_fct( v(mV), cai(mM)) {
:    h_inf = 1 / (1+exp((v+68.9)/6.5))
    h_inf = 1 / (1+exp((v+88.9+vshift)/7.8)) : Rodrigues and Oertel, VCN TStellate, 2006
    tau_s = 0.25*exp((v+183.6+vshift)/15.24) / tadj
    tau_f = 0.25*exp((v+158.6+vshift)/11.2) / (1+exp((v+75+vshift)/5.5)) / tadj

    alpha_s = h_inf / tau_s
    alpha_f = h_inf / tau_f
    beta_s = (1-h_inf) / tau_s
    beta_f = (1-h_inf) / tau_f

    C = cai*cai/(cac*cac)
    k2 = 4e-4 * tadj
    }
UNITSON
