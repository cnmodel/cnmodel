#!/usr/bin/env python

import numpy as np
from optparse import OptionParser
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui

class Table():
    def __init__(self):
        self.F = 0.4 # release probability (constant; no facilitiation)
        self.k0 = 1/1.75 # /s, baseline recovery rate from depletion (slow rate)
        self.kmax = 1/0.025 # /s, maximal recovery rate from depletion (fast rate)
        self.F_flag = 1.0 # compute facilitation
        self.dense = False  # include desensitization ? (not appropriate for some models)
        self.td = 0.05 #  time constant for calcium-dependent recovery
        self.kd =  0.7 # affinity of fast recovery proces for calcium sensor

        self.ts = 0.015 # decay time constant of glutatme clearance
        self.ks = 1000 # affinity of receptor desensitization for glutatmate
        # The large value means no desense occurs (otherwise, ks should be about
        # 0.6)

        self.kf = 0.6 # affinity of facilitiation process
        self.tf = 0.01 # make facilitation VERY slow

        self.dD = 0.02 # sets Ca that drives recovery(Ca influx per AP)
        # 0.02 yields rate-dep recovery in 100-300 Hz
        self.dF = 0.02 # sets Ca that drives facilitation

        self.glu = 0.3
        return

    def celltype(self, celltype):
        """
        Return the table that corresponds to the particular cell type
        """
        if celltype == 'bushy_epsc':
            table = self.bushy_epsc()
        elif celltype == 'bushy_ipsc':
            table = self.bushy_ipsc()
        elif celltype == 'stellate_epsc':
            table = self.stellate_epsc()
        elif celltype == 'stellate_ipsc':
            table = self.stellate_ipsc()
        elif celltype == 'XuF_epsc':
            table = self.XuF_epsc()
        elif celltype == 'default':
            table = self
        else:
            raise ValueError('Table: Celltype %s not implemented', celltype)
        return table

    def bushy_epsc(self):
        """ data is average of 3 cells studied with recovery curves and individually fit """
        self.F = 0.29366
        self.F_flag = 1.0 # compute facilitation
        self.desense = False  # include desensitization ? (not appropriate for some models)
        self.k0 = 0.52313
        self.kmax = 19.33805
        self.td = 0.01516
        self.kd = 0.11283
        self.ts = 17.9122
        self.ks = 11.531
        self.kf = 17.78
        self.tf = 0.00975
        self.dD = 0.57771
        self.dF = 0.60364
        self.glu = 2.12827
        return(self)

    def XuF_epsc(self):
        """ data is from XuFriedman  """
        self.F = 0.3
        self.F_flag = .00  # did not compute facilitation
        self.desense = True  # include desensitization ? (not appropriate for some models)
        self.k0 = 0.45
        self.kmax = 18
        self.td = 0.035
        self.kd = 0.07
        self.ts = 0.015
        self.ks = 0.6
        self.kf = 17.78
        self.tf = 0.00975
        self.dD = 0.57771
        self.dF = 0.60364
        self.glu = 1.0
        return(self)

    def stellate_epsc(self):
        """ data is average of 2 cells studied with recovery curves and individually fit """
        self.F = 0.30514
        self.desense = False  # include desensitization ? (not appropriate for some models)
        self.F_flag = 1.0 # compute facilitation
        self.k0 = 0.04045
        self.kmax = 38.4408
        self.td = 0.00549
        self.kd = 0.11911
        self.ts = 23.626
        self.ks = 21.06690
        self.kf = 12.24439
        self.tf = 0.01256
        self.dD = 2.07174
        self.dF = 2.08471
        self.glu = 8.298
        return(self)

    def bushy_ipsc(self):
        """ data is average of 4 cells studied with recovery curves and individually fit """
        self.F = 0.23382
        self.F_flag = 1.0 # compute facilitation
        self.desense = False  # include desensitization ? (not appropriate for some models)
        self.k0 = 0.67554
        self.kmax = 52.93832
        self.td = 0.08195
        self.kd = 0.28734
        self.ts = 0.17500
        self.ks = 4.57098
        self.kf = 16.21564
        self.tf = 0.12336
        self.dD = 2.21580
        self.dF = 1.17146
        self.glu = 1.90428
        return(self)

    def stellate_ipsc(self):
        """ data is average of 3 cells studied with recovery curves and individually fit """
        self.F = 0.23047
        self.F_flag = 1.0 # compute facilitation
        self.desense = False  # include desensitization ? (not appropriate for some models)
        self.k0 = 1.23636
        self.kmax = 45.34474
        self.td = 0.09809
        self.kd = 0.01183
        self.ts = 17.61450
        self.ks = 17.88618
        self.kf = 19.11424
        self.tf = 0.03228
        self.dD = 2.52072
        self.dF = 2.33317
        self.glu = 3.06948
        return(self)

class DKR():

    def __init__(self, callmode=None, table=None, stim=None, recovery=None, plot=False, celltype='bushy_epsc'):
        """
        # implementation eqs 1-4 of xu-f 2008///
        # or you could say, of Dittman et al, 2000
        # or Dittman and Regehr, 1998....
        #
        # call:
        # DKR no arguments just runs with some default parameters
        # DKR (mode, varargs):
        # if mode = 1, the table is passed in the second argument, and the plot is
        # generated.
        # if mode = 2, the table is in the second, the train in the third and the
        # recovery times in the 4th. For modeling/fitting
        # if mode = 3, the table is in the second, the train in the third and we
        # assume that there is no recovery data.
        # P. Manis 6/2008.
        # converted to Python 11/2009, 7/2014
        #
        """
        xout = []
        yout = []
        plotsvn = False
        self.plot = plot
        xtr = []
        if callmode == None:
            theTable = Table()
            table = theTable.celltype(celltype)
            #print table.kf, table.F
            freqlist = [50, 100, 200, 333, 400] # in Hz
            modelist = 0 # if 0, use regular; otherwise use exponential.
            plotvsn = True # control plots
            traindur = 0.5 # sec
            recovery = np.array([10.0,20.,30.,40.,50.,100.,200.,300.,500.,1000.]) # recovery points
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
            self.freq_run(mode, freqlist, traindur, stim, xtr, recovery, table)

    def freq_run(self, mode, freqlist, traindur, stim, xtr, recovery, table):
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
            pulse_train, rec = self.make_stimpattern(mode, freq, traindur, stim)
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
            S[0] = 1.0
        E[0] = Fn[0] * D[0]*S[0]  # initial EPSC variable
        # now for each stimulus pulse/release event
        for i in range(1, n_ptt):
            dt = pulse_train[i]-pulse_train[i-1]  # get time since last release
            [D[i], CaDi[i], Fn[i], CaFi[i]] = self.dstep(dt, table, D[i-1], CaDi[i-1], Fn[i-1], CaFi[i-1])
            glu[i] = (glu[i-1] + table.F * D[i] )* np.exp(-dt/table.ts)  # compute updated glutamate
            if table.desense:
                S[i] = table.ks/(table.ks + glu[i])
            else:
                S[i] = 1.0
            E[i] = Fn[i] * D[i] * S[i]  # Total release and desensitization.
        return D, S, E, Fn, glu, CaDi, CaFi

    def dstep(self, dt, T, Di, CaDi, Fi, CaFi):
        """
        [CaDi'; CaFi']'
        This is the discreet solver for each time step for the analytical equations as
        described in the papers.
        This version, including F, is from Dittman et al. 2000.
        calculate next D from Equations 15/16, Dittman et al. 2000
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

    def make_stimpattern(self, mode, freq, traindur, stim, n=0, xtr=None):
        rec = 0
        pt = None
        if mode == 0:
            pt = (1.0/freq)*np.ones(freq*traindur) # pulse train, msec.
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
    plot = True
    DKR(plot=plot, celltype = model)

