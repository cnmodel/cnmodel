: 6-state bushy cell model from Yang, Ramamurthy, Neef and Xu-Friedman
: "Low somatic sodium conducdtance... " J. Neurosci 26: 11999, 2016
:
: Implementation by P.B. Manis
: Kinetic data taken from Table 1
: Transition states from Figure 3A

NEURON {
    THREADSAFE
    SUFFIX nabu
    USEION na READ ena WRITE ina
    RANGE gna, gbar, vshift
    RANGE alpha, beta, gamma, gammao, deltao
    RANGE alpha_V, beta_V, gamma_V, gammao_V, deltao_V
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
}

: initialize parameters

PARAMETER {
    gbar = 0     (S/cm2)

    : alpha = 97.7196 (ms)
    alpha = 0.01023336 (/ms)
    alpha_V = 30.8128 (mV)
    
   : beta = 1.08888 (ms)
    beta = 0.9183748 (/ms)
    beta_V = -8.99497 (mV)
    
  :  gamma = 1.71473 (ms)
    gamma = 0.58318219 (/ms)
    gamma_V = 15.0341 (mV)
    
   : delta = 0.00164334 (ms)
    delta = 6085.1680 (/ms)
    delta_V = -13.7795 (mV)
    
   : gammao = 3.32014 (ms)
    gammao = 0.30119211 (/ms)
    gammao_V = 110.382 (mV)
    
   : deltao = 0.104053 (ms)
    deltao = 9610.4869 (/ms)
    deltao_V = -11.7251 (mV)
      
    vshift = 0 (mV)
}
    


ASSIGNED {
    v    (mV)
    ena  (mV)
    gna   (S/cm2)
    ina  (milliamp/cm2)
    alp   (ms)
    bet  (ms)
    gam   (ms)
    del  (ms)
    delo   (ms)
    gamo (ms)
}

STATE { C2 C1 O I2 I1 IO}

BREAKPOINT {
    SOLVE kfunc METHOD sparse
    gna = gbar * O
    ina = (gna)*(v - ena)
}

INITIAL { SOLVE kfunc METHOD sparse }

KINETIC kfunc {
    rates(v)
    ~ C1 <-> C2 (bet, 2*alp)
    ~ C1 <-> O (alp, 2*bet)
    
    ~ I1 <-> I2 (bet, 2*alp)
    ~ I1 <-> IO (alp, 2*bet)
    
    ~ I2 <-> C2 (del, gam)
    ~ I1 <-> C1 (del, gam)
    ~ IO <-> O (delo, gamo)
    
    CONSERVE C1 + C2 + I1 + I2 + IO + O = 1 
}


PROCEDURE rates(v(millivolt)) {
    LOCAL vS, maxrate
    maxrate = 8.00e+03 (/ms) : limiting value for reaction rates
                                     : See Patlak, 1991
    vS = v + vshift
	
    alp = (1./alpha)*exp(vS/alpha_V)
    bet = (1./beta)*exp(vS/beta_V)

    gam = (1./gamma)*exp(vS/gamma_V)
    del = (1./delta)*exp(vS/delta_V)

    
    delo = (1./deltao)*exp(vS/deltao_V)
    gamo = (1./gammao)*exp(vS/gammao_V)

	:alp = alp*maxrate / (alp+maxrate)
	:bet = bet*maxrate / (bet+maxrate)
	:gam = gam*maxrate / (gam+maxrate)
	:del = del*maxrate / (del+maxrate)
	:gamo = gamo*maxrate / (gamo+maxrate)
	:delo = delo*maxrate / (delo+maxrate)
}