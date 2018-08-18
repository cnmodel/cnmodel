#!/usr/bin/python


"""
Test adex model

"""

import numpy as np
from neuron import h
import neuron
import matplotlib.pyplot as plt
from collections import OrderedDict
import weave
from .gif.Filter_Rect_LogSpaced import *

class AdEx():
    def __init__(self):
        pass
        self.dt = 0.025
        self.tstop = 250.
        self.eta     = Filter_Rect_LogSpaced()    # nA, spike-triggered current (must be instance of class Filter)
        self.gamma   = Filter_Rect_LogSpaced()    # mV, spike-triggered movement of the firing threshold (must be instance of class Filter)
        
    def runone(self, pars, dt=0.025, tstop=250):
    
        self.cell = self.create_Adex_model(pars)
        self.dt = 0.025
        self.tstop = tstop
        
    #    for inj in np.linspace(pars['I']*0.5, pars['I']*2.0, 3):
        inj = pars['I']
        stim = h.Vector()
        self.Vm = h.Vector()
        self.Vm.record(self.cell(0.5)._ref_Vm_AdEx, sec=self.cell)
        istim = h.iStim(0.5, sec=self.cell)
        istim.delay = 5.
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0
        stim = np.zeros(int(tstop/h.dt))
        stim[int(20/h.dt):int(220/h.dt)] = inj
        cmd = h.Vector(stim)
        #cmd.play(istim._ref_i, h.dt, 0, sec=cell)
        cmd.play(self.cell(0.5)._ref_is_AdEx, h.dt, 0, sec=self.cell)

        rtime = h.Vector()
        rtime.record(h._ref_t)
        h.finitialize()
        h.t = 0.
        h.dt = self.dt
        while h.t < self.tstop:
            h.fadvance()

        tb = np.array(rtime)
        vcell = self.Vm.to_python()
        
        return(tb, vcell)

    def run_gifAI(self, pars, dt=0.025, tstep=250.):
        # Model parameters
        self.simulate(pars['I'], -65., pars)
        
        
        
    def plot_trace(self, i, panel, tb, vcell):
        axi = self.ax[i]
        axi.plot(tb, vcell)
        axi.set_ylim([-100, 10])
        axi.set_title(panel)  
        

    def create_Adex_model(self, pars):
        cell = h.Section()
        cell.L = 50

        cell.insert('AdEx')

        cell(0.5).cm_AdEx = pars['cm']
        cell(0.5).gl_AdEx = pars['gl']
        cell(0.5).el_AdEx = pars['el']
        cell(0.5).vt_AdEx = pars['Vt']
        cell(0.5).delt_AdEx = pars['dt']
        cell(0.5).a_AdEx = pars['a']
        cell(0.5).tauw_AdEx = pars['tauw']
        cell(0.5).b_AdEx = pars['b']
        cell(0.5).vr_AdEx = pars['Vr']
        cell(0.5).refract_AdEx=0.025
        return cell

    def create_GIF_model(self, pars):
        cell = h.Section()
        cell.L = 50
        cell.insert('GIF')
        cell(0.5).cm_GIF = pars['cm']
        cell(0.5).gl_GIF = pars['gl']
        cell(0.5).el_GIF = pars['el']
        cell(0.5).vt_GIF = pars['Vt_star']
        cell(0.5).vr_GIF = pars['Vr']
        cell(0.5).refract_GIF=pars['Tref']
        cell(0.5).DV_GIF = pars['DV']
        cell(0.5).lambda0_GIF = pars['lambda0']
        cell(0.5).b_GIF = pars['b']
        return cell        
        
        
    #
    # The following are from Gif.py...
    def simulate(self, I, V0, pars):
 
        """
        Simulate the spiking response of the GIF model to an input current I (nA) with time step dt.
        V0 indicate the initial condition V(0)=V0.
        The function returns:
        - time     : ms, support for V, eta_sum, V_T, spks
        - V        : mV, membrane potential
        - eta_sum  : nA, adaptation current
        - V_T      : mV, firing threshold
        - spks     : ms, list of spike times 
        """
 
        # Input parameters
        p_T         = len(I)
        p_dt        = self.dt
        
        # Model parameters
        p_gl        = pars['gl']
        p_C         = pars['cm'] 
        p_El        = pars['el']
        p_Vr        = pars['Vr']
        p_Tref      = pars['Tref']
        p_Vt_star   = pars['Vt_star']
        p_DV        = pars['DV']
        p_lambda0   = pars['lambda0']
        
        # Model kernels   
        (p_eta_support, p_eta) = self.eta.getInterpolatedFilter(self.dt)   
        p_eta       = p_eta.astype('double')
        p_eta_l     = len(p_eta)

        (p_gamma_support, p_gamma) = self.gamma.getInterpolatedFilter(self.dt)   
        p_gamma     = p_gamma.astype('double')
        p_gamma_l   = len(p_gamma)
      
        # Define arrays
        V = np.array(np.zeros(p_T), dtype="double")
        I = np.array(I, dtype="double")
        spks = np.array(np.zeros(p_T), dtype="double")                      
        eta_sum = np.array(np.zeros(p_T + 2*p_eta_l), dtype="double")
        gamma_sum = np.array(np.zeros(p_T + 2*p_gamma_l), dtype="double")            
 
        # Set initial condition
        V[0] = V0
         
        code =  """
                #include <math.h>
                
                int   T_ind      = int(p_T);                
                float dt         = float(p_dt); 
                
                float gl         = float(p_gl);
                float C          = float(p_C);
                float El         = float(p_El);
                float Vr         = float(p_Vr);
                int   Tref_ind   = int(float(p_Tref)/dt);
                float Vt_star    = float(p_Vt_star);
                float DeltaV     = float(p_DV);
                float lambda0    = float(p_lambda0);
           
                int eta_l        = int(p_eta_l);
                int gamma_l      = int(p_gamma_l);
                
                                                  
                float rand_max  = float(RAND_MAX); 
                float p_dontspike = 0.0 ;
                float lambda = 0.0 ;            
                float r = 0.0;
                
                                                
                for (int t=0; t<T_ind-1; t++) {
    
    
                    // INTEGRATE VOLTAGE
                    V[t+1] = V[t] + dt/C*( -gl*(V[t] - El) + I[t] - eta_sum[t] );
               
               
                    // COMPUTE PROBABILITY OF EMITTING ACTION POTENTIAL
                    lambda = lambda0*exp( (V[t+1]-Vt_star-gamma_sum[t])/DeltaV );
                    p_dontspike = exp(-lambda*(dt/1000.0)); // since lambda0 is in Hz, dt must also be in 1/Hz (this is why dt/1000.0)
                          
                          
                    // PRODUCE SPIKE STOCHASTICALLY
                    r = rand()/rand_max;
                    if (r > p_dontspike) {
                                        
                        if (t+1 < T_ind-1)                
                            spks[t+1] = 1.0; 
                        
                        t = t + Tref_ind;    
                        
                        if (t+1 < T_ind-1) 
                            V[t+1] = Vr;
                        
                        
                        // UPDATE ADAPTATION PROCESSES     
                        for(int j=0; j<eta_l; j++) 
                            eta_sum[t+1+j] += p_eta[j]; 
                        
                        for(int j=0; j<gamma_l; j++) 
                            gamma_sum[t+1+j] += p_gamma[j] ;  
                        
                    }
               
                }
                
                """
 
        vars =['p_T','p_dt','p_gl','p_C','p_El','p_Vr',
                'p_Tref','p_Vt_star','p_DV','p_lambda0',
                'V','I','p_eta','p_eta_l','eta_sum',
                'p_gamma','gamma_sum','p_gamma_l','spks' ]
        
        v = weave.inline(code, vars)

        time = np.arange(p_T)*self.dt
        
        eta_sum   = eta_sum[:p_T]     
        V_T = gamma_sum[:p_T] + p_Vt_star
     
        spks = (np.where(spks==1)[0])*self.dt
    
        return (time, V, eta_sum, V_T, spks)


#Type C (pF) gL (nS) EL (mV) VT (mV) dT (mV) a (nS) tauw (ms) b (pA) Vr (mV) I (pA)
# Values in the I column are hand-picked to try to reproduce the traces shown in Naud et al.,
# 2008. The values in Ims are the ones from the manuscript. 
        
    def fromfigure(self):
        
        tab=r"""Ty cm gl el   Vt dt a tauw b    Vr  I Ims
        4a        200 10 -70 -50 2 2   30  0   -58 300 500
        4b 200 12 -70 -50 2 2   300 60  -58 360 500
        4c 130 18 -58 -50 2 4   150 120 -50 250 400
        4d 200 10 -58 -50 2 2   120 100 -46 150 210
        4e 200 12 -70 -50 2 -10 300 0   -58 120 300
        4f 100 10 -65 -50 2 -10 90  30  -47  90 350
        4g 100 20 -70 -45 2 50   20  0   -64 900 110
        4h 100 12 -60 -50 2 -11 130 30  -48 120 160
        """
        lx = tab.splitlines()
        keys = lx[0].split()
        d = OrderedDict()
        for i in range(1,len(lx)):
            dline = lx[i].split()
            print (dline)
            if len(dline) == 0:
                continue
            fk = dline[0]
            d[fk] = {}
            for j, k in enumerate(keys[1:]):
                d[fk][k] = float(dline[j+1])
    
        f, self.ax = plt.subplots(len(list(d.keys()))/2, 2)
        self.ax = self.ax.ravel()
        for i, fig in enumerate(d.keys()):
            (tb, v) = self.runone(d[fig], self.ax[i])
            self.plot_trace(i, fig, tb, v)

#Type C (pF) gL (nS) EL (mV) VT (mV) dT (mV) a (nS) tauw (ms) b (pA) Vr (mV) I (pA)
# Values in the I column are hand-picked to try to reproduce the traces shown in Naud et al.,
# 2008. The values in Ims are the ones from the manuscript. 
# -------------------------
# GIF model parameters:
# -------------------------
# tau_m (ms):    21.237
# R (MOhm):    128.630
# C (nF):        0.165
# gl (nS):    0.007774
# El (mV):    -58.987
# Tref (ms):    4.000
# Vr (mV):    -32.586
# Vt* (mV):    -45.367
# DV (mV):    1.158
# -------------------------
#
    def showGIF(self):

        f, self.ax = plt.subplots(1, 3)
        self.ax = self.ax.ravel()
        pars = {
            'gl': 7.77,  # nS
            'cm': 165.,  # nF
            'el': -58.98,  # mV
            'Vr': -70,  # mV
            'Tref': 4.0,  # ms
            'Vt_star': -45.3,  # mV
            'DV': 1.158,
            'lambda0': 1.0,
            'Istim': 1000., # pA
            'dt': self.dt, # msec integration step
            }
        tstop = 250.
        npts = int(self.tstop/self.dt)
        pars['I'] = np.zeros(npts)
        pars['I'][int(npts*0.25):int(0.75*npts)] = pars['Istim']
        (time, V, eta_sum, V_T, spks) = self.simulate(pars['I'], pars['el'], pars)
        self.plot_trace(0, 'GIFAI', time, V)
        print('Array lengths: time = {0:d}, etasum = {1:d}, V_T = {2:d}'.format(len(time), len(eta_sum), len(V_T)))
        self.ax[1].plot(time[:len(eta_sum)], eta_sum, 'b')
        self.ax[2].plot(time[:len(V_T)], V_T, 'c')

if __name__ == '__main__':
    M = AdEx()
    #M.fromfigure()
    
    M.showGIF()
    
    plt.show()