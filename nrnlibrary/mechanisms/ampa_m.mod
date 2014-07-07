TITLE Model of AMPA receptors

COMMENT
-----------------------------------------------------------------------------

    Kinetic model of AMPA receptors
    ===============================

    6-state gating model:
    (scheme 1 from Raman and Trussell, Neuron 9:173-186, 1992)
        2 open states provide dual exponential response.
  
               O1
               |
    C -- C1 -- C2 -- O2
               |
               D

-----------------------------------------------------------------------------

  This mod file does not include mechanisms for the release and time course
  of transmitter; it is to be used in conjunction with a sepearate mechanism
  to describe the release of transmitter and that provides the concentration
  of transmitter in the synaptic cleft (to be connected to pointer C here).

  Default parameters are set for a miniature EPSC.

  Code based on Destexhe's ampa5.mod

B. Graham, Dept. of Computing Science & Maths, University of Stirling
(Contact: b.graham@cs.stir.ac.uk)
(previously IANC, Division of Informatics, University of Edinburgh)

CNS 2000 Version (19/11/02)


--------
Paul Manis (Otolaryngology/HNS and Cell and Molecular Physiology,
UNC Chapel Hill. contact: pmanis@med.unc.edu)
:
3/15/2005 Modifications:
1. Added Q10/qfac to allow temperature scaling. All rates in the state model are 
changed by the same factor. A Q10 of 1.5 gives a decay tau (single exponential fit
using Praxis algorithm in NEURON; using ampa_kinetics.hoc) of about 850 usec at 22
deg C and 570 usec at 33 deg C. These are consistent with the Raman and Trussell 1992
measurements in avians. The 850 usec is a bit fast for an EPSC, and could probably
be tuned by adjustment of some of the parameters below.

2. Brought several variables out to global (rather than range) so that we can
change them - Q10 and gmax in particular. note that gmax is in pS.

Max open probability is less than unity, so a gmax of 2500 yields 100 pA at
-60 mV. Therefore scaling by mini size must take this into account.

3/28/2005 Paul B. Manis
Added rectification to AMPA R.
Rectification is controlled by polyamine-style block of receptor. See
Donevan and Rogawski, 1995; Washburn et al., 1997. The equations used here
are from Washburn et al. 
The values given in the equation at the break point were determined from EPSCs in
5 21-d old DBA mice. Blocker = 45 (uM), Kd = 31.32, zd = 1.029.
Note that this should also reduce the maximal conductance.
Mode: if 1, use rectifying; if 0, use non-rectifying. Default is 1

-----------------------------------------------------------------------------
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    POINT_PROCESS AMPA_M
    POINTER CGLU
    RANGE C0, C1, C2, D, O1, O2
    RANGE g, rb,  Erev, gmax
    GLOBAL Rb, Ru1, Ru2, Rd, Rr, Ro1, Rc1, Ro2, Rc2
    GLOBAL vmin, vmax
	RANGE Q10, Mode
	RANGE zd, Kd0, PA
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (umho) = (micromho)
    (mM) = (milli/liter)
    (uM) = (micro/liter)
}

PARAMETER {

    Erev    = 7    (mV) : reversal potential
    gmax    = 500  (pS)  : maximal conductance
    vmin = -120 (mV)
    vmax = 100  (mV)
	Q10 = 1.5 : temperature sensitivity
    Mode = 1 : flag to control rectification calculation
: polyamine block parameters
	zd = 1.032
	PA = 45
	Kd0 = 31.e-6
	
: Rates

    Rb  = 13   (/mM /ms): binding 
                : diffusion limited (DO NOT ADJUST)
    Ru1 = 0.3  (/ms)    : unbinding (1st site)
    Ru2 = 200  (/ms)    : unbinding (2nd site)      
    Rd  = 30.0   (/ms)  : desensitization (WAS30.0)
    Rr  = 0.02 (/ms)    : resensitization 
    Ro1 = 100    (/ms)  : opening (fast)
    Rc1 = 2    (/ms)  : closing
    Ro2 = 2    (/ms)   : opening (slow)
    Rc2 = 0.25    (/ms) : closing
	aflag = 1 : Flag for control of printout of initial values.....
	
	
}

ASSIGNED {
    v       (mV)    : postsynaptic voltage
    i       (nA)    : current = g*(v - Erev)
    g       (pS)    : conductance
:	g0		(pS)	: conductance for voltage-dependent block by polyamines
	gvdep  (pS)   : voltage-dependence of conductance 
    CGLU       (mM)    : pointer to glutamate concentration
    rb      (/ms)   : binding
	qfac   : q10 factor for rate scaling

}

STATE {
    : Channel states (all fractions)
    C0      : unbound
    C1      : single glu bound
    C2      : double glu bound
    D       : single glu bound, desensitized
    O1      : open state 1
    O2      : open state 2
}

INITIAL {
    usetable = 0
	C0=1
    C1=0
    C2=0
    D=0
    O1=0
    O2=0
	qfac = Q10^((celsius-22)/10)
:	if(aflag == 1) {
:	VERBATIM
:	fprintf(stdout, "AMPA.MOD!!! gmax: %f    Q10 = %f  celsius = %f\n", gmax, Q10, celsius);
:ENDVERBATIM
:aflag = 0
:}
gvdepcalc(v)
	}

BREAKPOINT {
    SOLVE kstates METHOD sparse
    gvdepcalc(v)
    g = gmax * (O1 + O2)
    if ( Mode == 1) {
:		g0 = 1.0 + 0.6*exp((v-50)/40)  : eq. 5 of Washburn et al., 1997, slightly modified
:		gv = g0*(1/(1+PA/(Kd0*exp(-zd*v/25.3))))
		i = (1e-6) * g * gvdep * (v - Erev)
	}
	else {
		i = (1e-6)*g*(v-Erev)
	}
}

KINETIC kstates {

	
    rb = Rb * CGLU 

    ~ C0  <-> C1  (rb*qfac,Ru1*qfac)
    ~ C1 <-> C2 (rb*qfac,Ru2*qfac)
    ~ C2 <-> D  (Rd*qfac,Rr*qfac)
    ~ C2 <-> O1 (Ro1*qfac,Rc1*qfac)
    ~ C2 <-> O2 (Ro2*qfac,Rc2*qfac)

    CONSERVE C0+C1+C2+D+O1+O2 = 1
}
LOCAL g0
PROCEDURE gvdepcalc(v) {
	TABLE gvdep DEPEND PA, Kd0, zd FROM -100 TO 100 WITH 200
	
		g0 = 1.0 + 0.6*exp((v-50)/40)  : eq. 5 of Washburn et al., 1997, slightly modified
		gvdep = g0*(1/(1+PA/(Kd0*exp(-zd*v/25.3))))

}
