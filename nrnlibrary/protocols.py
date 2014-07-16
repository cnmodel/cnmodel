from neuron import h
#from neuron import *
#from nrnlibrary.pynrnutilities import *
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
        
    def run(self, ivrange, cell, durs=None, sites=None, reppulse=None, temp=22,
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
        istim.delay = 5.
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0
        (secmd, maxt, tstims) = make_pulse(stim)
        tend = maxt


        # Calculate current pulse levels
        iv_nstepi = int(np.ceil((imax - imin) / istep))
        iv_mini = imin
        iv_maxi = imax
        istep = (iv_maxi - iv_mini) / iv_nstepi
        iv_nstepi = iv_nstepi + 1
        for i in range(iv_nstepi):
            icur.append(float(i * istep) + iv_mini)
        nsteps = iv_nstepi
        
        self.current_cmd = np.array(icur)
        self.dt = dt
        self.temp = temp
        vec = {}

        for i in range(nsteps):
            # Generate current command for this level 
            stim['amp'] = icur[i]
            (secmd, maxt, tstims) = make_pulse(stim)
            vec['i_stim'] = h.Vector(secmd)
            
            # Connect recording vectors
            self['v_soma'] = cell.soma(0.5)._ref_v
            self['i_inj'] = istim._ref_i
            self['time'] = h._ref_t
            # connect current command vector
            vec['i_stim'].play(istim._ref_i, h.dt, 0, sec=cell.soma)

            # GO
            #h('secondorder=0')  # direct call fails; let hoc do the work
            h.dt = dt
            h.celsius = temp
            h.tstop = tend
            cell.cell_initialize()  # initialize the cell to it's rmp
            if i == 0:
                Rin, tau, v = cell.measure_rintau(auto_initialize=False)
                print '    *** Rin: %9.0f  tau: %9.1f   v: %6.1f' % (Rin, tau, v)
            h.frecord_init()
            #print('After initialization at {0:4.1f}: \n  vm0: {1:9.5f}  soma.v: {2:9.5f}'.format(
            #    temp, cell.vm0, cell.soma(0.5).v))
            #h.run()  # this seems to reset vm to -65 mV (a default), need to use fadvance
            while h.t < h.tstop:
                h.fadvance()

            k1 = int(stim['delay']/h.dt - 1)
#            print ('   V before step: {0:9.5f}'.format(self['v_soma'][k1]))
            self.voltage_traces.append(self['v_soma'])
            self.current_traces.append(self['i_inj'])
            self.time_values = np.array(self['time'])

    def peak_vm(self, window=0.5):
        """
        :param window: fraction of trace to look at to find peak value
        Return peak membrane voltage for each trace.

        """
        Vm = self.voltage_traces
        Icmd = self.current_cmd
        steps = len(Icmd)
        peakStart = self.durs[0] / self.dt
        peakStop = peakStart + (self.durs[1]*window) / self.dt # peak can be in first half
        Vpeak = []
        for i in range(steps):
            if Icmd[i] > 0:
                Vpeak.append(Vm[i][peakStart:peakStop].max())
            else:
                Vpeak.append(Vm[i][peakStart:peakStop].min())
        return np.array(Vpeak)
    
    def steady_vm(self, window=0.1):
        """
        :param window: fraction of window to use for steady-state measurement, taken
        immediately before the end of the step
        Return steady-state membrane voltage for each trace.
        """
        Vm = self.voltage_traces
        steps = len(Vm)
        steadyStop = (self.durs[0] + self.durs[1]) / self.dt
        steadyStart = steadyStop - (self.durs[1]*window) / self.dt  # measure last 10% of trace
        Vsteady = [Vm[i][steadyStart:steadyStop].mean() for i in range(steps)]
        return np.array(Vsteady)

    def spike_times(self, threshold=-40):
        """
        Return an array of spike times for each trace.
        """
        Vm = self.voltage_traces
        steps = len(Vm)
        spikes = []
        for i in range(steps):
            #dvdt = np.diff(Vm[i]) / self.dt
            #mask = (dvdt > 40).astype(int)
            mask = (Vm[i] > threshold).astype(int)
            indexes = np.argwhere(np.diff(mask) == 1)[:, 0] + 1
            times = indexes.astype(float) * self.dt
            spikes.append(times)
        return spikes

    def spike_filter(self, spikes, window=(0., np.inf)):
        """
        filter the spikes to only those occurring in a defined window.
        Required to compute input resistance in traces with no spikes during
        the stimulus, because some traces will have anodal break spikes.
        :param spikes: the list of spike trains returned from the spike_times method
        :param window: the window over which to look for spikes (in msec: default is
        the entire trace).

        return the spikes in a list
        """
        filteredspikes = []
        for i in range(len(spikes)):
            winspikes = []  # spikes is arranged by current; so this is for one level
            for j in range(len(spikes[i])):
                if spikes[i][j] >= window[0] and spikes[i][j] <= window[1]:
                    winspikes.append(spikes[i][j])
            filteredspikes.append(winspikes)  # now build filtered spike list
        return filteredspikes

    def rest_vm(self):
        """
        Return the resting membrane potential.
        """
        d = int(self.durs[0] / self.dt)
        return self.voltage_traces[-1][d//2:d].mean()
    
    def input_resistance(self, vmin=-70, imax=0):
        """
        Estimate resting input resistance.
        :param vmin: minimum voltage to use in computation
        :param imax: maximum current to use in computation.
        Return (slope, intercept) of linear regression for subthreshold traces
        near rest.
        Include traces in which spikes appear only AFTER the pulse using spike_filter
        """
        Vss = self.steady_vm()
        Icmd = self.current_cmd
        rawspikes = self.spike_times()
        spikes = self.spike_filter(rawspikes, window=[0., self.durs[0]+self.durs[1]])
        steps = len(Icmd)
        
        nSpikes = np.array([len(s) for s in spikes])
        # find traces with Icmd < 0, Vm > -70, and no spikes.
        vmask = Vss >= vmin
        imask = Icmd <= imax
        smask = nSpikes > 0
        mask = vmask & imask & ~smask
        if mask.sum() < 2:
            print('{0:<15s}: {1:s}'.format('vss', ', '.join(['{:.2f}'.format(v) for v in Vss])))
            print('{0:<15s}: {1:s}'.format('vmask', repr(vmask.astype(int))))
            print('{0:<15s}: {1:s} '.format('imask', repr(imask.astype(int))))
            print('{0:<15s}: {1:s}'.format('spikemask', repr(smask.astype(int))))
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
        print "\nMembrane resistance (chord): %0.1f MOhm" % s
        ivals = np.array([Icmd.min(), Icmd.max()])
        vvals = s * ivals + i
        line = pg.QtGui.QGraphicsLineItem(ivals[0], vvals[0], ivals[1], vvals[1])
        line.setPen(pg.mkPen(255, 0, 0, 70))
        line.setZValue(-10)
        IVplot.addItem(line, ignoreBounds=True)
        
        print "Resting membrane potential: %0.1f mV\n" % self.rest_vm()
        
        self.win = win

    
class VCCurve(Protocol):
    def __init__(self):
        super(VCCurve, self).__init__()

    def reset(self):
        super(VCCurve, self).reset()
        self.voltage_traces = []
        self.current_traces = []
        self.durs = None  # durations of current steps
        self.voltage_cmd = None # Current command levels
        self.time_values = None
        self.dt = None

    def run(self, vcrange, cell):
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
        self.reset()
        try:
            (vmin, vmax, vstep) = vcrange  # unpack the tuple...
        except:
            raise TypeError("run_iv argument 1 must be a tuple (imin, imax, istep)")

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
        vstep = (iv_maxv - iv_minv) / iv_nstepv
        for i in range(iv_nstepv):
            vcmd.append(float(i * vstep) + iv_minv)
        nreps = iv_nstepv
        for i in range(nreps):
            # Connect recording vectors
            self['v_soma'] = cell.soma(0.5)._ref_v
            self['i_inj'] = istim._ref_i
            self['time'] = h._ref_t
            vstim.amp2 = vcmd[i]
            h.tstop = tend
            h.init()
            h.run()
            self.voltage_traces.append(self['v_soma'])
            self.current_traces.append(self['i_inj'])
            self.time_values = np.array(self['time'])


    def steady_im(self, window=0.1):
        """
        :param window: fraction of window to use for steady-state measurement, taken
        immediately before the end of the step
        Return steady-state membrane current for each trace.
        """
        Im = self.current_traces
        steps = len(Im)
        steadyStop = (self.durs[0] + self.durs[1]) / self.dt
        steadyStart = steadyStop - (self.durs[1]*window) / self.dt
        Isteady = [Im[i][steadyStart:steadyStop].mean() for i in range(steps)]
        return np.array(Isteady)


# def run_democlamp(cell, dend, vsteps=[-60,-70,-60], tsteps=[10,50,100]):
#     """
#     Does some stuff.
#
#     """
#     f1 = pylab.figure(1)
#     gs = GS.GridSpec(2, 2,
#                        width_ratios=[1, 1],
#                        height_ratios=[1, 1])
#
#     # note numbering for insets goes from 1 (upper right) to 4 (lower right)
#     # counterclockwise
#     pA = f1.add_subplot(gs[0])
#     pAi = INSETS.inset_axes(pA, width="66%", height="40%", loc=2)
#     pB = f1.add_subplot(gs[1])
#     pBi = INSETS.inset_axes(pB, width="66%", height="40%", loc=4)
#     pC = f1.add_subplot(gs[2])
#     pCi = INSETS.inset_axes(pC, width="66%", height="40%", loc=2)
#     pD = f1.add_subplot(gs[3])
#     pDi = INSETS.inset_axes(pD, width="66%", height="40%", loc=1)
#     #h.topology()
#
#     Ld = 0.5
#     Ld2 = 1.0
#
#     VClamp = h.SEClamp(0.5, cell)
#     VClamp.dur1 = tsteps[0]
#     VClamp.amp1 = vsteps[0]
#     VClamp.dur2 = tsteps[1]
#     VClamp.amp2 = vsteps[1]
#     VClamp.dur3 = tsteps[2]
#     VClamp.amp3 = vsteps[2]
#     Rs0 = 10.
#     VClamp.rs = Rs0
#     compensation = [0., 70., 95.]
#     cms = [cell.cm*(100.-c)/100. for c in compensation]
#
#     vrec = h.iStim(Ld, sec=dend[0])
#     vrec.delay = 0
#     vrec.dur = 1e9 # these actually do not matter...
#     vrec.iMax = 0.0
#     vrec2 = h.iStim(Ld2, sec=dend[0])
#     vrec2.delay = 0
#     vrec2.dur = 1e9 # these actually do not matter...
#     vrec2.iMax = 0.0
#
#     stim = {}
#     stim['NP'] = 1
#     stim['Sfreq'] = 20 # stimulus frequency
#     stim['delay'] = tsteps[0]
#     stim['dur'] = tsteps[1]
#     stim['amp'] = vsteps[1]
#     stim['PT'] = 0.0
#     stim['hold'] = vsteps[0]
# #    (secmd, maxt, tstims) = make_pulse(stim)
#     tend = 79.5
#     linetype = ['-', '-', '-']
#     linethick = [0.5, 0.75, 1.25]
#     linecolor = [[0.66, 0.66, 0.66], [0.4, 0.4, 0.3], 'k']
#     n = 0
#     vcmds = [-70, -20]
#     vplots = [(pA, pAi, pC, pCi), (pB, pBi, pD, pDi)]
#     for m,  VX in enumerate(vcmds):
#         stim['amp'] = VX
#         pl = vplots[m]
#         print m, VX
#         (secmd, maxt, tstims) = make_pulse(stim)
#         for n, rsc in enumerate(compensation):
#             vec={}
#             for var in ['VCmd', 'i_inj', 'time', 'Vsoma', 'Vdend',
#                         'Vdend2', 'VCmdR']:
#                 vec[var] = h.Vector()
#             VClamp.rs = Rs0*(100.-rsc)/100.
#             cell.cm = cms[n]
#            # print VClamp.rs, cell.cm, VClamp.rs*cell.cm
#             vec['VCmd'] = h.Vector(secmd)
#             vec['Vsoma'].record(cell(0.5)._ref_v, sec=cell)
#             vec['Vdend'].record(dend[0](Ld)._ref_v, sec=dend[0])
#             vec['time'].record(h._ref_t)
#             vec['i_inj'].record(VClamp._ref_i, sec=cell)
#
#             vec['VCmdR'].record(VClamp._ref_vc, sec=cell)
#             VClamp.amp2 = VX
#             #            vec['VCmd'].play(VClamp.amp2, h.dt, 0, sec=cell)
#
#             h.tstop = tend
#             h.init()
#             h.finitialize(-60)
#             h.run()
#             vc = np.asarray(vec['Vsoma'])
#             tc = np.asarray(vec['time'])
#
#             # now plot the data, raw and as insets
#             for k in [0, 1]:
#                 pl[k].plot(vec['time'], vec['i_inj'], color=linecolor[n], linestyle = linetype[n], linewidth=linethick[n])
#                 yl = pl[k].get_ylim()
#                 if k == 0:
#                     pass
#                     #pl[k].set_ylim([1.5*yl[0], -1.5*yl[1]])
#                 else:
#                     pass
#
#             for k in [2,3]:
#                 pl[k].plot(vec['time'], vec['Vsoma'], color=linecolor[n], linestyle = linetype[n], linewidth=linethick[n])
#                 pl[k].plot(vec['time'], vec['VCmdR'], color=linecolor[n], linestyle = '--', linewidth=1, dashes=(1,1))
#                 pl[k].plot(vec['time'], vec['Vdend'], color=linecolor[n], linestyle = linetype[n], linewidth=linethick[n], dashes=(3,3))
#                 if VX < vsteps[0]:
#                     pl[k].set_ylim([-72, -40])
#                 else:
#                     pl[k].set_ylim([-62,VX+30])
#
#     ptx = 10.8
#     pBi.set_xlim([9.8, ptx])
#     pAi.set_xlim([9.8, ptx])
#     PH.setX(pAi, pCi)
#     PH.setX(pBi, pDi)
#     pD.set_ylim([-65, 10])
# #    PH.setY(pC, pCi) # match Y limits
#     PH.cleanAxes([pA, pAi, pB, pBi, pC, pCi, pD, pDi])
#     PH.formatTicks([pA, pB, pC, pD], axis='x', fmt='%d')
#     PH.formatTicks([pC, pD], axis='y', fmt='%d')
#     PH.calbar(pAi, [ptx-1, 0, 0.2, 2.])
#     PH.calbar(pCi, [ptx-1, -50., 0.2, 10])
#     PH.calbar(pBi, [ptx-1, 0, 0.2, 10])
#     PH.calbar(pDi, [ptx-1, -50., 0.2, 20])
#     pylab.draw()
#     pylab.show()
