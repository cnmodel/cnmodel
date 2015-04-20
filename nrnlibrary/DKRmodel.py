#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DKRmodel computes synaptic release dynamics for arbitrary input trains using the model
 of Dittman, J.S., Krieger, A.C. and Regehr, W.G.  "Interplay between Facilitation,
 Depression, and Residual Calcium at Three Presynaptic Terminals",

 The Journal of Neuroscience, February 15, 2000, 20(4):1374–1385

 The model is modified to include the receptor desensitization term from

 H Yang and MA Xu-Friedman, Relative roles of different mechanisms of
 depression at the mouse endbulb of Held.
 J Neurophysiol 2008, 99: 2510-2521.
 The desensitization term is controlled by a flag and normally should be False unless attempting
 to replicate the Yang and Xu-Friedman results. (Desensitization is a function of receptors, and
 should be part of the receptor model.)

There are 2 classes in this file:
    Class Table sets the kinetic parameters of the model. A variety of parameters are provided
    corresonding to both average and individual measurements from bushy and planar multipolar
    cells of the CBA mouse VCN, collected by Ruili Xie in 2008-2009, measured in 2.5 mM Ca,
    1.5 mM Mg (see Xie and Manis, J. Neurosci. "Target-Specific IPSC Kinetics Promote Temporal
    Processing in Auditory Parallel Pathways" The Journal of Neuroscience, January 23, 2013 • 33(4):1598 –1614
    These are from simultaneous fits of the kinetic model to a range of stimulus
    frequencies from 50 to 400 Hz.

    Class DKR solves the differential equations based on the event times. It assumes a stimulus
    train of 20 pulses, with recovery measurements (should be generalized).

TODO:
    1. generalize DKR to accept arbitrary input stimulus train (partially implemented, not tested)
    2. Allow better access of stimulus parameters for regular trains by user.
    3. Incorporate fitting algorithm (new class here) to fit raw data to model and return parameters.
        This is already done in the original matlab code, just need to incorporate here.

Author: Paul B. Manis, Ph.D., UNC Chapel Hill


"""
import numpy as np
from optparse import OptionParser
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui

class Table():
    def __init__(self):
        self.dkr_constant = False  # If True, kf / dF is constant
        self.F_flag = 1.0  # compute facilitation (if 0, F is set to 1; see Yang and Xu-Friedman)
        self.desense = False  # include desensitization (True to emulate Yang and Xu-Friedman)
                            # (generally should be in incorporated into a receptor model, but may be useful
                            # for fits to systems with desensitizing receptors)
        self.F = None  # F1 in dkr. resting release probability (constant) after long period, no stim
        self.rho = None  # Not in original model
        self.k0 = None  # /s, baseline recovery rate from depletion (slow rate)
        self.kmax = None  # /s, maximal recovery rate from depletion (fast rate)
        self.td = None  #  time constant for calcium-dependent recovery
        self.kd = None  # affinity of fast recovery process for calcium sensor

        self.ts = None  # decay time constant of glutamate clearance
        self.ks =  None  # affinity of receptor desensitization for glutamate
        # The large value means no desensitization occurs (otherwise, ks should be about
        # 0.6)

        self.kf = None  # affinity of facilitation process
        self.tf = None  # make facilitation VERY slow

        self.dD = None  # sets Ca that drives recovery(Ca influx per AP)
        # 0.02 yields rate-dep recovery in 100-300 Hz
        self.dF = None  # sets Ca that drives facilitation

        self.glu = None  # glutamate concentration in cleft (estimated from model)
        return

    def celltype(self, celltype, select='average'):
        """
        Return the table that corresponds to the particular cell type and
        selected cell or average
        """
        if celltype == 'bushy_epsc':
            self.bushy_epsc(select)
        elif celltype == 'bushy_ipsc':
            self.bushy_ipsc(select)
        elif celltype == 'stellate_epsc':
            self.stellate_epsc(select)
        elif celltype == 'stellate_ipsc':
            self.stellate_ipsc(select)
        elif celltype == 'XuF_epsc':
            self.XuF_epsc()
        elif celltype == 'DKRSC':
            self.DKR_SC()
        elif celltype == 'default':
            pass
        else:
            raise ValueError('Table: Celltype %s not implemented', celltype)
        # print check for consistentcy in table. In X&M2013, we fitted dD and dF,
        # but they are actually constrained as follows,
        # "Kf/dF is entirely determined by experimentally observed value of facilitation
        #  and the initial release probability"
        # (restate Eq 7 DKR00):

        return(self)

    def bushy_epsc(self, select = "average"):
        """ data is from Xie and Manis, 2013: average of 3 cells studied with recovery curves and individually fit """
        # mean 3 cells: 12sep08b, 19sep08a, 02sep08b
        # Synfit_Results_EPSCs_IPSCs.xlsx

        self.F_flag = 1.0 # compute facilitation
        self.desense = False  # include desensitization ? (not appropriate for some models)
        if select == 'average':
            self.F = 0.29366
            self.k0 = 0.52313
            self.kmax = 19.33805
            self.kd = 0.11283
            self.ks = 11.531
            self.kf = 17.78
            self.td = 0.01516
            self.ts = 17.912
            self.tf = 0.00975
            self.dD = 0.57771
            self.dF = 0.60364
            self.glu = 2.12827

        elif select == 'nodepression':  # hand tuned model
            print 'no depression'
            self.F = 0.2
            self.k0 = 0.45 # 0.52313
            self.kmax = 18 # 19.33805
            self.kd = 0.7 # 0.01 # 0.11283
            self.ks = 0.6 # 11.531
            self.kf = 5. # 17.78
            self.ts = 0.015 # 17.912
            self.tf = 0.1 # 0.00975
            self.td = self.tf/2. # 0.01516 #  0.01516
            self.dD = 0.5
            self.dF = 0.60364
            self.glu = 2.12827

        # self.F = 0.29366
        # self.k0 = 0.52313  #/ 1000.
        # self.kmax = 19.33805  #/ 1000.
        # self.kd = 0.11283
        # self.ks = 11.531
        # self.kf = 17.78
        # self.td = 0.01516
        # self.ts = 17.9122
        # self.tf = 0.00975
        # self.dD = 0.57771
        # self.dF = 0.60364
        # self.glu = 2.12827
        #
        elif select == '12Sep08b':
            self.F = 0.3177
            self.k0 = 1.09962
            self.kmax = 19.29
            self.kd = 0.04849
            self.ks = 6.08
            self.kf = 10.13
            self.td = 0.01396
            self.ts = 12.27
            self.tf = 0.01125
            self.dD = 0.49382
            self.dF = 1.38095
            self.glu = 0.903

        elif select in ['19Sep08a', 'exemplar']:
            # # 19Sep08a  (Exemplar in figure 8 of Xie and Manis, 2013).
            self.F = 0.26686
            self.k0 = 0.22446
            self.kmax = 18.87
            self.kd = 0.15201
            self.ks = 14.92727
            self.kf = 30.0
            self.td = 0.01456
            self.ts = 29.32
            self.tf = 0.00900
            self.dD = 1.04385
            self.dF = 0.41597
            self.glu = 5.4805

        elif select in ['02Sep08b']:
            # # 19Sep08a  (Exemplar in figure 8 of Xie and Manis, 2013).
            self.F = 0.29544
            self.k0 = 0.24530
            self.kmax = 19.85092
            self.kd = 0.13801
            self.ks = 13.58577
            self.kf = 13.20655
            self.td = 0.01456
            self.ts = 12.14891
            self.tf = 0.00900
            self.dD = 0.19547
            self.dF = 0.0140
            self.glu = 0.001

        else:
            raise ValueError ("Bushy cell selection %s not known", select)

    def stellate_epsc(self, select='average'):
        """ data is from Xie and Manis, 2013 """
#        self.F = 0.434

        self.desense = False  # include desensitization ? (not appropriate for some models)
        self.F_flag = 1.0  # compute facilitation

        if select == 'average':
            # Mean, 3 cells (15sep08c, 08sep08c 29dec08b)
            # from Synfit_Results_EPSCs_IPSCs.xlsx
            self.F = 0.43435
            self.k0 = 0.0672
            self.kmax = 52.82713
            self.td = 0.00398
            self.kd = 0.08209
            self.ts = 16.917
            self.ks = 14.24460
            self.kf = 18.16292
            self.tf = 0.01138
            self.dD = 2.46535
            self.dF = 1.44543
            self.glu = 5.86564

        elif select == '15Sep08c':
            # from Synfit_Results_EPSCs_IPSCs.xlsx
            self.F = 0.17371
            self.k0 = 0.03636
            self.kmax = 28.4114
            self.td = 0.00983
            self.kd = 0.23733
            self.ts = 19.50927
            self.ks = 11.39048
            self.kf = 20.57235
            self.tf = 0.01613
            self.dD = 0.57484
            self.dF = 1.16943
            self.glu = 0.001

        elif select == '08Sep08c':
            # from Synfit_Results_EPSCs_IPSCs.xlsx
            self.F = 0.43658
            self.k0 = 0.04454
            self.kmax = 48.46876
            self.td = 0.00115
            self.kd = 0.00090
            self.ts = 27.74209
            self.ks = 30.74333
            self.kf = 3.91642
            self.tf = 0.00900
            self.dD = 3.56865
            self.dF = 3.0000
            self.glu = 16.59592

        elif select == '29Dec08c':
            # from Synfit_Results_EPSCs_IPSCs.xlsx
            self.F = 0.69276
            self.k0 = 0.12060
            self.kmax = 81.60123
            self.td = 0.00096
            self.kd = 0.00806
            self.ts = 3.50000
            self.ks = 0.60000
            self.kf = 30.0
            self.tf = 0.00900
            self.dD = 3.25256
            self.dF = 0.16658
            self.glu = 1.0

        else:
            raise ValueError ("Stellate cell EPSC selection %s not known", select)

    def bushy_ipsc(self, select='average'):
        """
        Paramegers are from fits to data, Xie and Manis, 2013
        """
        self.F_flag = 1.0 # compute facilitation
        self.desense = False  # include desensitization ? (not appropriate for some models)
        if select == 'average':
            self.F = 0.18521
            self.k0 = 2.29700
            self.kmax = 27.6667
            self.td = 0.12366
            self.kd = 0.12272
            self.ts = 9.59624
            self.ks = 8.854469
            self.kf = 5.70771
            self.tf = 0.37752
            self.dD = 4.00335
            self.dF = 0.72605
            self.glu = 5.61985
        # self.k0 = 2.30
        # self.kmax = 27.7
        # self.td = 0.00124
        # self.kd = 0.123
        # self.ts = 0.00960
        # self.ks = 8.85
        # self.kf = 5.71
        # self.tf = 0.00378
        # self.dD = 4.00
        # self.dF = 0.726
        # self.glu = 5.62
        if select == '30Aug08f':
            #30aug08f: facilitates at low F noted bad fit
            self.F_flag = 1.0 # compute facilitation
            self.desense = False  # include desensitization ? (not appropriate for some models)
            self.F = 0.19575
            self.k0 = 0.03636
            self.kmax = 21.4477
            self.td = 0.00248
            self.kd = 0.00090
            self.ts = 0.175
            self.ks = 5.57537
            self.kf = 0.34234
            self.tf = 0.15151
            self.dD = 3.70
            self.dF = 0.04413
            self.glu = 20.0

        elif select == '30Aug08b':
            #30aug08b: noted good fit
            self.F_flag = 1.0 # compute facilitation
            self.desense = False  # include desensitization ? (not appropriate for some models)
            self.F = 0.28734
            self.k0 = 0.44069
            self.kmax = 14.88
            self.td = 0.00263
            self.kd = 0.00090
            self.ts = 30.68
            self.ks = 38.54
            self.kf = 1.756
            self.tf = 0.32375
            self.dD = 3.70
            self.dF = 0.13624
            self.glu = 19.254

        else:
            raise ValueError ("Bushy cell IPSC selection %s not known", select)

    def stellate_ipsc(self, select='average'):
        """ data is average of 3 cells studied with recovery curves and individually fit """
        # self.F = 0.23047
        self.F_flag = 1.0 # compute facilitation
        self.desense = False  # include desensitization ? (not appropriate for some models)
        # self.k0 = 1.23636
        # self.kmax = 45.34474
        # self.td = 0.09809
        # self.kd = 0.01183
        # self.ts = 17.61450
        # self.ks = 17.88618
        # self.kf = 19.11424
        # self.tf = 0.03228
        # self.dD = 2.52072
        # self.dF = 2.33317
        # self.glu = 3.06948
        if select == 'average':
            # mean 3 cells from Synfit_Results_EPSCs_IPSCs.xlsx
            self.F = 0.23047
            self.k0 = 1.23636 #/ 1000.0
            self.kmax = 45.34474 #/ 1000.0
            self.td = 98.09
            self.kd = 0.01183
            self.ts = 17614.50
            self.ks = 17.88618
            self.kf = 19.11424
            self.tf = 32.28
            self.dD = 2.52072
            self.dF = 2.33317
            self.glu = 3.06948

        else:
            raise ValueError ("Stellate cell IPSC selection %s not known", select)

    # Test cases of other published data:
    def XuF_epsc(self):
        """ Parameters from from Yang and XuFriedman  """
        self.F_flag = .0  # did not compute facilitation
        self.desense = True  # include desensitization ? (not appropriate for some models)
        self.F = 0.3
        self.k0 = 0.45
        self.kmax = 18.
        self.td = 0.035
        self.kd = 0.07
        self.ts = 0.015
        self.ks = 0.6
        self.kf = 17.78
        self.tf = 0.00975
        self.dD = 0.57771
        self.dF = 0.60364
        self.glu = 1.0

    def DKR_SC(self):
        """ Parameters from DKR Schaffer collateral (Table 2)  """
        self.dkr_constant = True
        self.F_flag = 1.0  #
        self.desense = False  # include desensitization ? (not appropriate for some models)
        self.F = 0.24
        self.rho = 2.2
        self.k0 = 2.0
        self.kmax = 30.0
        self.td = 0.050  # in seconds, not ms
        self.kd = 2.0
       # self.ts = 0.015  # not used
       # self.ks = 0.6  # not used
       # self.kf = 17.78
        self.tf = 0.100  # in seconds, not ms
        self.dF = 0.5
        #rho = (1.0 - self.F)*self.F2/self.F
        #self.dD = 0.5 # they never give the values in the paper...
        self.glu = 1.0   # arbitrary here


class DKR():

    def __init__(self, callmode=None, table=None, stim=None, recovery=None,
                 plot=False, celltype='bushy_epsc', select='average'):
        """
        implementation eqs of Dittman et al, 2000
        with modification by Yang and Xu-Friedman for desensitization parameter.

        call:
        DKR no arguments just runs with some default parameters
        DKR (mode, varargs):
        if mode = 1, the table is passed in the second argument, and the plot is
        generated.
        if mode = 2, the table is in the second, the train in the third and the
        recovery times in the 4th. For modeling/fitting
        if mode = 3, the table is in the second, the train in the third and we
        assume that there is no recovery data.
        P. Manis 6/2008.
        converted to Python 11/2009, 7/2014
        """
        xout = []
        yout = []
        plotsvn = False
        self.plot = plot
        xtr = []
        if callmode == None:
            theTable = Table()
            table = theTable.celltype(celltype, select=select)
            print table
            freqlist = [50, 100, 200, 400] # in Hz
            modelist = 0 # if 0, use regular; otherwise use exponential.
            plotvsn = True # control plots
            traindur = None # 0.5 # sec
            trainpulses = 20  # number of pulses in train
            recovery = np.array([10.0,20.,30.,40.,50.,100.,200.,300.,500.,1000., 1500., 2000.]) # recovery points
            recovery = recovery*0.001 # convert to seconds
        else:
            if callmode in (0,1):
                freqlist = [50,100,200,400] # in Hz
                if stim is not None:
                    modelist = stim
                else:
                    modelist = 0 # if 0, use regular; otherwise use exponential.
                if callmode is 1:
                    plotvsn = True # control plots

                traindur = 0.5 # sec
                recovery = 0.001*np.array([10.0,20,30,40,50,100,200,300,500,1000]) # recovery points
            elif callmode == 2:# with recovery data
                xtr = recovery # and recoveries (matched)
                freqlist = np.zeros(len(stim))
                for k in range(0, len(stim)):
                    freqlist[k] = 1.0/np.mean(np.diff(stim[k]))
                modelist = 2
            elif callmode == 3: # no recovery
                xtr = [] # and recoveries (matched)
                freqlist = np.arange((0,len(stim))*100)
                modelist = 2
            else:
                raise ValueError('DKR - call mode not recognized\n')
        mode = modelist
        # if mode == 1:
        #     ntrial = 10
        #     for trial in range(len(ntrial)):
        #         for n, fr in enumerate(freqlist):
        #             self.one_run(mode, freqlist[n], traindur, stim, xtr, recovery, table)
        #             if trial == 0:
        #                 self.ER = self.ESum[len(pulse_train):]
        #             else:
        #                 self.ER = self.ER + self.ESum[len(pulse_train):]
        #     self.ER = self.ER/ntrial  # averate over trials
        if mode != 1:
            ntrial = 1
            self.freq_run(mode, freqlist, traindur, trainpulses, stim, xtr, recovery, table)

    def freq_run(self, mode, freqlist, traindur, trainpulses, stim, xtr, recovery, table):
        """
        Compute the model responses for many frequencies with multiple recovery times
        """
        # lists of results, each array holding data for one frequency
        self.ESum = []  # EPSC summary, normalized
        self.DSum = []  # Depression
        self.SSum = []  # deSensitization
        self.CaDiSum = []  # calcium in depression
        self.CaFiSum = []  # calcium for faciitiation
        self.FnSum = []  # facilitation
        self.gluSum = []  # glutmate (concentration)
        self.PulseTrain = []  # pulse train times

        for n, freq in enumerate(freqlist):
            pulse_train, rec = self.make_stimpattern(mode, freq, traindur, trainpulses, stim)
            D, S, E, Fn, glu, CaDi, CaFi = self.compute_recovery(pulse_train, recovery, table)
            self.ESum.append(E) #  = np.append(self.ESum, E[-1])
            self.DSum.append(D) #  = np.append(self.DSum, D[-1])
            self.SSum.append(S)  #= np.append(self.SSum, S[-1])
            self.CaDiSum.append(CaDi)  # = np.append(self.CaDiSum, CaDi[-1])
            self.CaFiSum.append(CaFi)  # = np.append(self.CaFiSum, CaFi[-1])
            self.FnSum.append(Fn) #  = np.append(self.FnSum, Fn[-1])
            self.gluSum.append(glu)  # = np.append(self.gluSum, glu[-1])
            self.PulseTrain.append(np.append(pulse_train, np.amax(pulse_train) + recovery))

        if self.plot is True: # plot if requested
            self.show(mode, freqlist)
        return

    def compute_recovery(self, pulse_train, recovery, table):
        """
        Compute the recovery function for a given pulse train and list of recovery times
        Returns:
        states of internal variables, collapsed across recovery times
        """
        nr = len(recovery)
        D = []
        S = [0]*nr
        E = [0]*nr
        Fn = [0]*nr
        glu = [0]*nr
        CaDi = [0]*nr
        CaFi = [0]*nr
        for j in range(len(recovery)):  # for each recovery pulse, compute the full response.
            ptt = np.append(pulse_train, np.amax(pulse_train) + recovery[j]) # pick just one recovery pulse per train to prevent interactions
            Dj, Sj, Ej, Fnj, gluj, CaDij, CaFij = self.compute_DKR(ptt, table)
            if j == 0:
                D = Dj
                S = Sj
                E = Ej
                Fn = Fnj
                glu = gluj
                CaDi = CaDij
                CaFi = CaFij
            else:
                D = np.append(D, Dj[-1]) # just add the last pulse
                S = np.append(S, Sj[-1])
                E = np.append(E, Ej[-1])
                Fn = np.append(Fn, Fnj[-1])
                glu = np.append(glu, gluj[-1])
                CaDi = np.append(CaDi, CaDij[-1])
                CaFi = np.append(CaFi, CaFij[-1])

        return D, S, E, Fn, glu, CaDi, CaFi

    def compute_DKR(self, pulse_train, table):
        """
        Compute Dittman et al. 2000 equations for one stimulus train (pt)
        """
        n_ptt = len(pulse_train)
        D = np.zeros(n_ptt) # depletion (1-D)
        S = np.zeros(n_ptt) # desens (1-S)
        E = np.zeros(n_ptt) # EPSP
        glu = np.zeros(n_ptt)
        CaDi = np.zeros(n_ptt)
        CaFi = np.zeros(n_ptt)
        Fn = np.zeros(n_ptt)
        Fn[0] = table.F
        glu[0] = table.glu
        dt = 1 # set initially to  a very long time since the last pulse
        # compute initial state
        (D[0], CaDi[0], Fn[0], CaFi[0]) = self.dstep(dt, table, 1, 0.1,  table.F, 0.1)
        if table.desense:  # do we include desensitization in this calculation or not?
            S[0] = table.ks/(table.ks + glu[0])
        else:
            S[0] = 1.0  # assume no desensitization in measured response.
        E[0] = Fn[0] * D[0]*S[0]  # initial EPSC variable
        # now for each stimulus pulse/release event
        for i in range(1, n_ptt):
            dt = pulse_train[i]-pulse_train[i-1]  # get time since last release
            [D[i], CaDi[i], Fn[i], CaFi[i]] = self.dstep(dt, table, D[i-1], CaDi[i-1], Fn[i-1], CaFi[i-1])
            # handle difference between the 2 models:
            if table.desense:
                glu[i] = (glu[i-1] + table.F * D[i] )* np.exp(-dt/table.ts)  # compute updated glutamate
                S[i] = table.ks/(table.ks + glu[i])
            else:
                glu[i] = (glu[i-1] + table.F * D[i] )  # assumes clearance is very fast
                S[i] = 1.0
            E[i] = Fn[i] * D[i] * S[i]  # Total release and desensitization.
        return D, S, E, Fn, glu, CaDi, CaFi

    def dstep(self, dt, T, Di, CaDi, Fi, CaFi):
        """
        [CaDi'; CaFi']'
        This is the discreet solver for each time step for the analytical equations as
        described in the papers.
        This version, including F, is from Dittman et al. 2000.
        Calculate next D from Equations 15/16, Dittman et al. 2000
        Handles two cases: known dD and dF, or the case where these are constant.
        Parameters:
        dt is time since last step
        Di is D at last time step
        Fi is F at the last time step
        T is the table (parameters/constants).
        CaDi is the calcium term for depression
        CaFi is the calcium term for facilitation

        returns:
        D, CaD, F, and CaF for the step advance
        """

        if T.dkr_constant:
            # set dD and dF based on rho
            #rho = (1.0 - self.F)*self.F2/self.F  # need a value for rho...

            T.kf = T.dF/(((1-T.F)/((T.F/(1-T.F)*T.rho-T.F)))-1)
            CaFi = CaFi + T.dF
            CaDi = CaFi  # assume equal, then test it.
        else:
            CaDi = CaDi + T.dD
            CaFi = CaFi + T.dF
        CaDn = CaDi*np.exp(-dt/T.td) # change with next step
        CaFn = CaFi*np.exp(-dt/T.tf)
        # r = 1/dt # convert to hz (units are SECONDS).
        eta = (T.kd/CaDi + 1)/(T.kd/CaDi + np.exp(-dt/T.td))
        eta = np.power(eta,(-(T.kmax-T.k0)*T.td))
        Dn = 1-(1-(1-Fi)*Di)*np.exp(-T.k0*dt)*eta
        Fn = T.F + T.F_flag*(1-T.F)/(1+T.kf/CaFn)
        return Dn, CaDn, Fn, CaFn

    def make_stimpattern(self, mode, freq, traindur, trainpulses, stim, n=0, xtr=None):
        rec = 0
        pt = None
        if mode == 0:
            if traindur is not None:
                pt = (1.0/freq)*np.ones(freq*traindur) # pulse train, msec.
            elif trainpulses is not None:
                pt = (1.0/freq)*np.ones(trainpulses)
            else:
                raise ValueError('make-stimpattern: either triandur or train pulses must be specified')
            pt = np.cumsum(pt)
            pt = 0.010+pt-pt[0] # make starts all the same.
        elif mode == 1:
            pt=[]
            pt[0] = np.randexpo(freq)
            k = 1
            while np.max(np.cumsum(pt)) < traindur:
                pt[k]=np.randexpo(freq) ##ok<AGROW>
                k = k + 1
            pt=np.cumsum(pt) # that's the stimulus train.
            pt = 0.010+pt-pt[0] # make starts all the same.
        elif mode == 2 and stim is not None: # get from incoming data
            pt = stim[n]
            if xtr is not None:
                rec = xtr[n]
            else:
                rec = 0
        else:
            pass
        return pt, rec

    def show(self, mode, freqlist):
        colors = ['r','g','b', 'c', 'm', 'y']
        # if mode == 0:
        #     mark_face = 'w' # regular are open symbols
        # elif mode == 2:
        #     mark_face = 'k'
        # else:
        #     mark_face = colors[n] # poisson are filled

        app = pg.mkQApp()
        win = pg.GraphicsWindow('Dittman-K-R Model')
        win.resize(1000, 800)
        self.p1 = win.addPlot(labels={'left': 'E (normalized)', 'bottom': 'Time (ms)'})
        self.p1.setLogMode(x=True)
        win.nextRow()
        self.p2 = win.addPlot(labels={'left': 'D', 'bottom': 'Time (ms)'})
        win.nextRow()
        self.p3 = win.addPlot(labels={'left': 'Fn', 'bottom': 'Time (ms)'})
        win.nextRow()
        self.p4 = win.addPlot(labels={'left': 'CaD (d) CaF (x)', 'bottom': 'Time (ms)'})
        leg = self.p4.addLegend(size=(40, 100), offset=(60, 20))
        #print dir(leg)
        win.nextRow()
        self.p5 = win.addPlot(labels={'left': 'S', 'bottom': 'Time (ms)'})
    #p1.plot(pt, E[0:-1]/esum[0],  colors[n] + 's-', markerfacecolor =mark_face, markersize =1.5)
        symsize = 6
        for fr in range(len(freqlist)):
            xout = self.PulseTrain[fr]
   #         print self.ESum[fr]
            self.p1.plot(xout, self.ESum[fr], pen = pg.mkPen(colors[fr]), symbolBrush=pg.mkBrush(colors[fr]),
                         symbol='o', symbolSize=symsize)
            self.p2.plot(xout, self.DSum[fr], pen=pg.mkPen(colors[fr]), symbolBrush=pg.mkBrush(colors[fr]),
                         symbol='o', symbolSize=symsize)
            self.p3.plot(xout, self.FnSum[fr], pen=pg.mkPen(colors[fr]), symbolBrush=pg.mkBrush(colors[fr]),
                         symbol='o', symbolSize=symsize)
            self.p4.plot(xout, self.CaFiSum[fr], pen=pg.mkPen(colors[fr]), symbolBrush=pg.mkBrush(colors[fr]),
                         symbol='o', symbolSize=symsize)
            self.p4.plot(xout, self.CaDiSum[fr], pen=pg.mkPen(colors[fr]),
                         symbol='+', symbolBrush=pg.mkBrush(colors[fr]),
                         symbolSize=symsize, name=str(freqlist[fr]))
            self.p5.plot(xout, self.SSum[fr], pen=pg.mkPen(colors[fr]), symbolBrush=pg.mkBrush(colors[fr]),
                         symbol='d', symbolSize=symsize)
        win.show()
        QtGui.QApplication.instance().exec_()


def test():
    DKR(plot=True, celltype = 'bushy_epsc')

if __name__ == "__main__":

    parser=OptionParser()
    parser.add_option("-b", "--bushyepsc", dest="be", action='store_true',
                      help="test bushy_epsc", default=False)
    parser.add_option("-s", "--stellateepsc", dest="se", action='store_true',
                      help="test stellate_epsc", default=False)
    parser.add_option("-B", "--bushyipsc", dest="bi", action='store_true',
                      help="test bushy_ipsc", default=False)
    parser.add_option("-S", "--stellateipsc", dest="si", action='store_true',
                      help="test bushy_epsc", default=False)
    parser.add_option("-F", "--XuF", dest="xuf", action='store_true',
                      help="test XuF measurements", default=False)
    parser.add_option("--DKRSC", dest="dkrsc", action='store_true',
                      help="test DKR Schaffer Collateral measurements", default=False)
    parser.add_option("-c", "--select", action="store", type="string", dest="select",
                      help='select cells or average', default='average')

    (options, args) = parser.parse_args()
    model = 'bushy_epsc'
    if options.be:
        model = 'bushy_epsc'
    if options.bi:
        model = 'bushy_ipsc'
    if options.se:
        model = 'stellate_epsc'
    if options.si:
        model = 'stellate_ipsc'
    if options.xuf:
        model = 'XuF_epsc'
    if options.dkrsc:
        model = 'DKRSC'

    plot = True
    DKR(plot=plot, celltype=model, select=options.select)

