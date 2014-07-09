import os
import os.path
from neuron import h
#from neuron import *
#from nrnlibrary.pynrnutilities import *
import pylibrary.Utility as U
import pylibrary.PlotHelpers as PH
import numpy as np
import scipy
import scipy.integrate
import scipy.stats

try:
    import pyqtgraph as pg
    HAVE_PG = True
except ImportError:
    HAVE_PG = False

from .util.stim import make_pulse

#import matplotlib as MP # must call first... before pylag/pyplot or backends
#MP.use('Qt4Agg')

#import matplotlib.gridspec as GS
#import mpl_toolkits.axes_grid1.inset_locator as INSETS
#import mpl_toolkits.axes_grid1.anchored_artists as ANCHOR

#stdFont = 'Arial'
#import  matplotlib.pyplot as pylab
#pylab.rcParams['interactive'] = False
#pylab.rcParams['mathtext.default'] = 'sf'
## next setting allows pdf font to be readable in Adobe Illustrator
#pylab.rcParams['pdf.fonttype'] = 42
#pylab.rcParams['figure.facecolor'] = 'white'



class Protocol(object):
    """
    Base class providing common tools for running, analyzing, and displaying
    simulations.
    """
    def __init__(self):
        pass

    def reset(self):
        self._vectors = {}
    
    def __setitem__(self, name, variable):
        """
        Record *variable* during the next run.
        """
        vec = h.Vector()
        self._vectors[name] = vec
        vec.record(variable)
        
    def __getitem__(self, name):
        """
        Return a np array for previously recorded data given *name*.
        """
        return np.array(self._vectors[name])


class IVCurve(Protocol):
    def __init__(self):
        super(IVCurve, self).__init__()
    
    def reset(self):
        super(IVCurve, self).reset()
        self.voltage_traces = []
        self.durs = None  # durations of current steps
        self.current_cmd = None # Current command levels
        self.current_traces = []
        self.time_values = None
        self.dt = None
        
    def run(self, ivrange, cell, durs=None, sites=None, reppulse=None, temp=32,
            dt=0.025):
        """
        Run a current-clamp I/V curve on *cell*.
        
        Parameters:
        ivrange : tuple
            (min, max, step)
        cell : Cell
            The Cell instance to test.
        durs : tuple
            durations of (pre, pulse, post) regions of the command
        sites : list
            Sections to add recording electrodes
        reppulse : 
            stimulate with pulse train
        temp : 
            temperature of simulation (32)
        dt : 
            timestep of simulation (0.025)
        """
        self.reset()
        
        try:
            (imin, imax, istep) = ivrange # unpack the tuple...
        except:
            raise TypeError("run_iv argument 1 must be a tuple (imin, imax, istep)")
        
        
        # Configure IClamp
        if durs is None:
            durs = [10.0, 100.0, 50.0]
            
        self.durs = durs
        
        icur = []
        # set up stimulation with a pulse train
        if reppulse is None:
            stim = {
                'NP': 1,
                'delay': durs[0],
                'dur': durs[1],
                'amp': 1.0,
                'dt': dt,
                }
        else:
            stim = {
                'NP': 10,
                'Sfreq': 50.0,
                'delay': 10.0,
                'dur': 2,
                'amp': 1.0,
                'PT': 0.0,
                'dt': dt,
                }
        istim = h.iStim(0.5, sec=cell.soma)
        istim.delay = 0
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0
        (secmd, maxt, tstims) = make_pulse(stim)
        tend = maxt


        # Calculate current pulse levels
        iv_nstepi = int(np.ceil((imax - imin) / istep))
        iv_mini = imin
        iv_maxi = imax
        nreps = iv_nstepi
        istep = (iv_maxi - iv_mini) / iv_nstepi
        iv_nstepi = iv_nstepi + 1
        for i in range(iv_nstepi):
            icur.append(float(i * istep) + iv_mini)
        nsteps = iv_nstepi
        
        self.current_cmd = np.array(icur)
        self.dt = dt
        self.temp = temp
        vec = {}
        
        splist = np.zeros(nsteps)
        meanVss = np.zeros(nsteps)
        meanIss = np.zeros(nsteps)
        minVpk = np.zeros(nsteps)
        
        for i in range(nsteps):
            # Set up recording vectors
            #for var in ['v_soma', 'i_inj', 'time', 'm', 'h', 'ah', 'bh', 'am',
                        #'bm', 'gh', 'ik', 'ina', 'inat', 'i_stim']:
                #vec[var] = h.Vector()
            #if sites is not None:
                #for j in range(len(sites)):
                    #vec['v_meas_%d' % (j)] = h.Vector()
                    
            # Generate current command for this level 
            stim['amp'] = icur[i]
            (secmd, maxt, tstims) = make_pulse(stim)
            vec['i_stim'] = h.Vector(secmd)
                
            
            # Connect recording vectors
            self['v_soma'] = cell.soma(0.5)._ref_v
            #self['ik'] = cell.soma(0.5)._ref_ik
            #natFlag = False
            #try:
                #vec['inat'].record(cell.soma(0.5)._ref_inat, sec=cell.soma)
                #natFlag = True
            #except:
                #vec['ina'].record(cell.soma(0.5)._ref_ina, sec=cell.soma)
                #pass
            
            #if sites is not None:
                #for j in range(len(sites)):
                    #if sites[j] is not None:
                        #vec['v_meas_%d' % (j)].record(
                            #sites[j](0.5)._ref_v, sec=sites[j])
            self['i_inj'] = istim._ref_i
            #vec['gh'].record(cell.soma().ihvcn._ref_i, sec=cell.soma)
            self['time'] = h._ref_t
            
            # connect current command vector
            vec['i_stim'].play(istim._ref_i, h.dt, 0, sec=cell.soma)

            # GO
            h.dt = dt
            h.celsius = temp
            h.tstop = tend
            h.init()
            h.run()
            
            self.voltage_traces.append(self['v_soma'])
            self.current_traces.append(self['i_inj'])
            self.time_values = np.array(self['time'])

    def peak_vm(self):
        """
        Return peak membrane voltage for each trace.
        """
        Vm = self.voltage_traces
        Icmd = self.current_cmd
        steps = len(Icmd)
        peakStart = self.durs[0] / self.dt
        peakStop = peakStart + 10. / self.dt
        Vpeak = []
        for i in range(steps):
            if Icmd[i] > 0:
                Vpeak.append(Vm[i][peakStart:peakStop].max())
            else:
                Vpeak.append(Vm[i][peakStart:peakStop].min())
        return np.array(Vpeak)
    
    def steady_vm(self):
        """
        Return steady-state membrane voltage for each trace.
        """
        Vm = self.voltage_traces
        steps = len(Vm)
        steadyStop = (self.durs[0] + self.durs[1]) / self.dt
        steadyStart = steadyStop - 30. / self.dt
        Vsteady = [Vm[i][steadyStart:steadyStop].mean() for i in range(steps)]
        return np.array(Vsteady)

    def spike_times(self):
        """
        Return an array of spike times for each trace.
        """
        Vm = self.voltage_traces
        steps = len(Vm)
        spikes = []
        for i in range(steps):
            #dvdt = np.diff(Vm[i]) / self.dt
            #mask = (dvdt > 40).astype(int)
            mask = (Vm[i] > -40.).astype(int)
            indexes = np.argwhere(np.diff(mask) == 1)[:, 0] + 1
            times = indexes.astype(float) * self.dt
            spikes.append(times)
        return spikes

    def rest_vm(self):
        """
        Return the resting membrane potential.
        """
        d = int(self.durs[0] / self.dt)
        return self.voltage_traces[-1][d//2:d].mean()
    
    def input_resistance(self, vmin=-70, imax=0):
        """
        Estimate resting input resistance.
        Return (slope, intercept) of linear regression for subthreshold traces
        near rest.
        """
        Vss = self.steady_vm()
        Icmd = self.current_cmd
        spikes = self.spike_times()
        steps = len(Icmd)
        
        nSpikes = np.array([len(s) for s in spikes])
        
        # find traces with Icmd < 0, Vm > -70, and no spikes.
        vmask = Vss >= vmin
        imask = Icmd <= imax
        smask = nSpikes == 0
        mask = vmask & imask & smask
        if mask.sum() < 2:
            print vmask.astype(int)
            print imask.astype(int)
            print smask.astype(int)
            raise Exception("Not enough traces to do linear regression.")
        
        # Use these to measure input resistance by linear regression.
        reg = scipy.stats.linregress(Icmd[mask], Vss[mask])
        (slope, intercept, r, p, stderr) = reg
        
        return slope, intercept

    def show(self):
        """
        Plot results from run_iv()
        """
        if not HAVE_PG:
            raise Exception("Requires pyqtgraph")
        
        #
        # Generate figure with subplots
        #
        app = pg.mkQApp()
        win = pg.GraphicsWindow()
        win.resize(1000, 800)
        Vplot = win.addPlot(labels={'left': 'Vm (mV)', 'bottom': 'Time (ms)'})
        rightGrid = win.addLayout(rowspan=2)
        win.nextRow()
        Iplot = win.addPlot(labels={'left': 'Iinj (nA)', 'bottom': 'Time (ms)'})
        
        IVplot = rightGrid.addPlot(labels={'left': 'Vm (mV)', 'bottom': 'Icmd (nA)'})
        IVplot.showGrid(x=True, y=True)
        rightGrid.nextRow()
        spikePlot = rightGrid.addPlot(labels={'left': 'Iinj (nA)', 'bottom': 'Spike times (ms)'})
        rightGrid.nextRow()
        FIplot = rightGrid.addPlot(labels={'left': 'Spike count', 'bottom': 'Iinj (nA)'})
        
        win.ci.layout.setRowStretchFactor(0, 10)
        win.ci.layout.setRowStretchFactor(1, 5)

        #
        # Plot simulation and analysis results
        #
        Vm = self.voltage_traces
        Iinj = self.current_traces
        Icmd = self.current_cmd
        t = self.time_values
        steps = len(Icmd)

    
        # plot I, V traces
        colors = [(i, steps*3./2.) for i in range(steps)]
        for i in range(steps):
            Vplot.plot(t, Vm[i], pen=colors[i])
            Iplot.plot(t, Iinj[i], pen=colors[i])


        # I/V relationships
        IVplot.plot(Icmd, self.peak_vm(), symbol='o')
        IVplot.plot(Icmd, self.steady_vm(), symbol='s')


        # F/I relationship and raster plot
        spikes = self.spike_times()
        for i,times in enumerate(spikes):
            spikePlot.plot(x=times, y=[Icmd[i]]*len(times), pen=None, 
                           symbol='d', symbolBrush=colors[i])
        FIplot.plot(x=Icmd, y=[len(s) for s in spikes], symbol='o')
        
        
        # Print Rm, Vrest 
        (s, i) = self.input_resistance()
        print "Membrane resistance: %0.1f MOhm" % s
        ivals = np.array([Icmd.min(), Icmd.max()])
        vvals = s * ivals + i
        line = pg.QtGui.QGraphicsLineItem(ivals[0], vvals[0], ivals[1], vvals[1])
        line.setPen(pg.mkPen(255, 0, 0, 70))
        line.setZValue(-10)
        IVplot.addItem(line, ignoreBounds=True)
        
        print "Resting membrane potential: %0.1f mV" % self.rest_vm()
        
        self.win = win

    

def run_vc(vmin, vmax, vstep, cell):
    """
    Run voltage-clamp I/V curve.
    
    Parameters:
    vmin : float
        Minimum voltage step value
    vmax : 
        Maximum voltage step value
    vstep :
        Voltage difference between steps
    cell :
        The Cell instance to test.
    """
    vstim = h.SEClamp(0.5, cell.soma) # use our new iclamp method
    vstim.dur1 = 50.0
    vstim.amp1 = -60
    vstim.dur2 = 500.0
    vstim.amp2 = -60.0
    vstim.dur3 = 400
    vstim.amp3 = -60.0
    vstim.rs = 0.01
    cell.soma.cm = 0.001
    vcmd = []
    tend = 900.0
    iv_nstepv = int(np.ceil((vmax - vmin) / vstep))
    iv_minv = vmin
    iv_maxv = vmax
    nreps = iv_nstepv
    vstep = (iv_maxv - iv_minv) / iv_nstepv
    for i in range(iv_nstepv):
        vcmd.append(float(i * vstep) + iv_minv)
#    tend = 160
    nreps = iv_nstepv
    vec = {}
    f1 = pylab.figure(1)
    gs = GS.GridSpec(2, 2,
                       width_ratios=[3, 1],
                       height_ratios=[3, 1])

    p1 = f1.add_subplot(gs[0])
    p2 = f1.add_subplot(gs[1])
    p3 = f1.add_subplot(gs[2])
    p4 = f1.add_subplot(gs[3])
#    p1 = f1.add_subplot(2,1,1)
#    p2 = f1.add_subplot(2,1,2)
#    p3 = f1.add_subplot(3,1,3)
    meanVss = np.zeros(nreps)
    meanIss = np.zeros(nreps)
    for i in range(nreps):
        for var in ['v_soma', 'i_inj', 'time', 'm', 'h',
                    'ah', 'bh', 'am', 'bm']:
            vec[var] = h.Vector()
        vstim.amp2 = vcmd[i]
        h.tstop = tend
        vec['v_soma'].record(cell.soma(0.5)._ref_v, sec=cell.soma)
        vec['i_inj'].record(vstim._ref_i, sec=cell.soma)
        vec['time'].record(h._ref_t)
        h.init()
        h.run()
#        tvec = arange(0, h.tstop, h.dt)
        p3.plot(vec['time'], vec['v_soma'])
        p1.plot(vec['time'], vec['i_inj'])
        (meanVss[i], r1) = U.measure(
            'mean', vec['time'], vec['v_soma'], 500, 550)
        (meanIss[i], r2) = U.measure(
            'mean', vec['time'], vec['i_inj'], 500, 550)


    p1.set_xlim(0, tend)
    p3.set_xlim(0, tend)
    p2.plot(meanVss, meanIss, color='r', linestyle='-', marker='s')
    PH.cleanAxes([p1, p2, p3, p4])
    PH.calbar(p1, [600, -0.7, 100., 0.1])
    PH.calbar(p3, [600, -85., 100., 10])
    PH.noaxes(p4)
    pylab.draw()
    pylab.show()
    
def run_democlamp(cell, dend, vsteps=[-60,-70,-60], tsteps=[10,50,100]):
    """
    Does some stuff.
    
    """
    f1 = pylab.figure(1)
    gs = GS.GridSpec(2, 2,
                       width_ratios=[1, 1],
                       height_ratios=[1, 1])

    # note numbering for insets goes from 1 (upper right) to 4 (lower right)
    # counterclockwise
    pA = f1.add_subplot(gs[0])
    pAi = INSETS.inset_axes(pA, width="66%", height="40%", loc=2)
    pB = f1.add_subplot(gs[1])
    pBi = INSETS.inset_axes(pB, width="66%", height="40%", loc=4)
    pC = f1.add_subplot(gs[2])
    pCi = INSETS.inset_axes(pC, width="66%", height="40%", loc=2)
    pD = f1.add_subplot(gs[3])
    pDi = INSETS.inset_axes(pD, width="66%", height="40%", loc=1)
    #h.topology()
    
    Ld = 0.5
    Ld2 = 1.0
    
    VClamp = h.SEClamp(0.5, cell)
    VClamp.dur1 = tsteps[0]
    VClamp.amp1 = vsteps[0]
    VClamp.dur2 = tsteps[1]
    VClamp.amp2 = vsteps[1]
    VClamp.dur3 = tsteps[2]
    VClamp.amp3 = vsteps[2]
    Rs0 = 10.
    VClamp.rs = Rs0
    compensation = [0., 70., 95.]
    cms = [cell.cm*(100.-c)/100. for c in compensation]
    
    vrec = h.iStim(Ld, sec=dend[0])
    vrec.delay = 0
    vrec.dur = 1e9 # these actually do not matter...
    vrec.iMax = 0.0
    vrec2 = h.iStim(Ld2, sec=dend[0])
    vrec2.delay = 0
    vrec2.dur = 1e9 # these actually do not matter...
    vrec2.iMax = 0.0

    stim = {}
    stim['NP'] = 1
    stim['Sfreq'] = 20 # stimulus frequency
    stim['delay'] = tsteps[0]
    stim['dur'] = tsteps[1]
    stim['amp'] = vsteps[1]
    stim['PT'] = 0.0
    stim['hold'] = vsteps[0]
#    (secmd, maxt, tstims) = make_pulse(stim)
    tend = 79.5
    linetype = ['-', '-', '-']
    linethick = [0.5, 0.75, 1.25]
    linecolor = [[0.66, 0.66, 0.66], [0.4, 0.4, 0.3], 'k'] 
    n = 0
    vcmds = [-70, -20]
    vplots = [(pA, pAi, pC, pCi), (pB, pBi, pD, pDi)]
    for m,  VX in enumerate(vcmds):
        stim['amp'] = VX
        pl = vplots[m]
        print m, VX
        (secmd, maxt, tstims) = make_pulse(stim)
        for n, rsc in enumerate(compensation):
            vec={}
            for var in ['VCmd', 'i_inj', 'time', 'Vsoma', 'Vdend',
                        'Vdend2', 'VCmdR']:
                vec[var] = h.Vector()
            VClamp.rs = Rs0*(100.-rsc)/100.
            cell.cm = cms[n]
           # print VClamp.rs, cell.cm, VClamp.rs*cell.cm
            vec['VCmd'] = h.Vector(secmd)
            vec['Vsoma'].record(cell(0.5)._ref_v, sec=cell)
            vec['Vdend'].record(dend[0](Ld)._ref_v, sec=dend[0])
            vec['time'].record(h._ref_t)
            vec['i_inj'].record(VClamp._ref_i, sec=cell)

            vec['VCmdR'].record(VClamp._ref_vc, sec=cell)
            VClamp.amp2 = VX
            #            vec['VCmd'].play(VClamp.amp2, h.dt, 0, sec=cell)

            h.tstop = tend
            h.init()
            h.finitialize(-60)
            h.run()
            vc = np.asarray(vec['Vsoma'])
            tc = np.asarray(vec['time'])
            
            # now plot the data, raw and as insets
            for k in [0, 1]:
                pl[k].plot(vec['time'], vec['i_inj'], color=linecolor[n], linestyle = linetype[n], linewidth=linethick[n])
                yl = pl[k].get_ylim()
                if k == 0:
                    pass
                    #pl[k].set_ylim([1.5*yl[0], -1.5*yl[1]])
                else:
                    pass
            
            for k in [2,3]:
                pl[k].plot(vec['time'], vec['Vsoma'], color=linecolor[n], linestyle = linetype[n], linewidth=linethick[n])
                pl[k].plot(vec['time'], vec['VCmdR'], color=linecolor[n], linestyle = '--', linewidth=1, dashes=(1,1))
                pl[k].plot(vec['time'], vec['Vdend'], color=linecolor[n], linestyle = linetype[n], linewidth=linethick[n], dashes=(3,3))
                if VX < vsteps[0]:
                    pl[k].set_ylim([-72, -40])
                else:
                    pl[k].set_ylim([-62,VX+30])

    ptx = 10.8
    pBi.set_xlim([9.8, ptx])
    pAi.set_xlim([9.8, ptx])
    PH.setX(pAi, pCi)
    PH.setX(pBi, pDi)
    pD.set_ylim([-65, 10])
#    PH.setY(pC, pCi) # match Y limits
    PH.cleanAxes([pA, pAi, pB, pBi, pC, pCi, pD, pDi])
    PH.formatTicks([pA, pB, pC, pD], axis='x', fmt='%d')
    PH.formatTicks([pC, pD], axis='y', fmt='%d')
    PH.calbar(pAi, [ptx-1, 0, 0.2, 2.])
    PH.calbar(pCi, [ptx-1, -50., 0.2, 10])
    PH.calbar(pBi, [ptx-1, 0, 0.2, 10])
    PH.calbar(pDi, [ptx-1, -50., 0.2, 20])
    pylab.draw()
    pylab.show()    
