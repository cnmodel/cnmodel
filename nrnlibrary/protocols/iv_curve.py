from neuron import h
import numpy as np
import scipy
import scipy.integrate
import scipy.stats
import scipy.optimize
import nrnlibrary.util.expfitting as expfitting

try:
    import pyqtgraph as pg
    HAVE_PG = True
except ImportError:
    HAVE_PG = False

from ..util.stim import make_pulse
from .protocol import Protocol


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
        self.cell = cell
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

        iextend = []
        if self.durs[2] > 50.:
            iextend = np.ones(int((self.durs[2]-50)/stim['dt']))
            secmd = np.append(secmd, secmd[-1]*iextend)
        tend = maxt + len(iextend)*stim['dt']


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
                self.custom_init()
                Rin, tau, v = cell.measure_rintau(auto_initialize=False)
                print '    *** Rin: %9.0f  tau: %9.1f   v: %6.1f' % (Rin, tau, v)

            cell.cell_initialize()  # initialize the cell to it's rmp
            self.custom_init()
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
        peakStart = int(self.durs[0] / self.dt)
        peakStop = int(peakStart + (self.durs[1]*window) / self.dt) # peak can be in first half
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
        steadyStop = int((self.durs[0] + self.durs[1]) / self.dt)
        steadyStart = int(steadyStop - (self.durs[1]*window) / self.dt)  # measure last 10% of trace
        Vsteady = [Vm[i][steadyStart:steadyStop].mean() for i in range(steps)]
        return np.array(Vsteady)

    def spike_times(self, threshold=None):
        """
        Return an array of spike times for each trace.
        
        :param threshold: Optional threshold at which to detect spikes. By 
        default, this queries cell.spike_threshold.
        """
        if threshold is None:
            threshold = self.cell.spike_threshold
        
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
    
    def input_resistance_tau(self, vmin=-70, imax=0):
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
            return 0., 0.
            #raise Exception("Not enough traces to do linear regression.")
        
        # Use these to measure input resistance by linear regression.
        reg = scipy.stats.linregress(Icmd[mask], Vss[mask])
        (slope, intercept, r, p, stderr) = reg

        # also measure the tau in the same traces:
        window = 0.5  # only use first half of trace for fitting.
        peakStart = int(self.durs[0] / self.dt)
        peakStop = int(peakStart + (self.durs[1]*window) / self.dt) # peak can be in first half

        fits = []
        tx = self.time_values[peakStart:peakStop] - self.durs[0]
        nexp = 1
        if nexp == 2:
            fitter = expfitting.ExpFitting(nexp=2, 
                    initpars={'dc': vmin, 'a1': -10., 't1': 5., 'a2': 2., 'delta': 3.0},
                    bounds={'dc': (-120, 0.), 'a1': (-50, 0.), 't1': (0.1, 25.), 'a2': (-50., 50.), 'delta': (3., 50.)}
            )
        elif nexp == 1:
            fitter = expfitting.ExpFitting(nexp=1, 
                    initpars={'dc': vmin, 'a1': 2., 't1': 5.},
                    bounds={'dc': (-120, 0.), 'a1': (-50, 0.), 't1': (0.1, 50.)}
            )
        for i, m in enumerate(mask):
            if m and (self.rest_vm() - Vss[i]) > 1:
                print( ' i: {0:3d}  v: {1:6.2f}'.format(i, Vss[i]))
                # 10/9/2014: replaced with lmfit version - much more stable and has bounds control
                #fitParams, fitCovariances = scipy.optimize.curve_fit(self.expFunc,
                #                                                     tx, self.voltage_traces[i][peakStart:peakStop],
                #                                                     p0 = [vmin, 2., 2., 5., 15.], maxfev = 5000)
                
                fitParams = fitter.fit(tx, self.voltage_traces[i][peakStart:peakStop], fitter.fitpars)
                fits.append(fitParams)
        tau1 = np.mean([f['t1'] for f in fits])
        dc = np.mean([f['dc'] for f in fits])
        print ('mean dc:   {:7.3f} '.format(dc))
        print ('mean tau1: {:7.3f}'.format(tau1))
        print ('dc:      {:s}'.format(''.join(['{0:7.3f}  '.format(f['dc'].value) for f in fits])))
        print ('amp1:    {:s}'.format(''.join([ '{0:7.3f}  '.format(a['a1'].value) for a in fits])))
        print ('tau1:    {:s}'.format(''.join([ '{0:7.3f}  '.format(a['t1'].value) for a in fits])))
        if nexp == 2:
            tau2 = np.mean([f['delta'].value*f['t1'].value for f in fits])
            print ('mean tau2:  {:7.3f}'.format(tau2))
            print ('amp2:    {:s}'.format(''.join(['{:7.3f}  '.format(a['a2'].value) for a in fits])))
            print ('tau2:    {:s}'.format(''.join(['{:7.3f}  '.format(a['t1'].value*a['delta'].value) for a in fits])))
        else:
            tau2 = 0.
        return slope, intercept, [tau1, tau2]


    def expFunc(self, t, d, a1, r1, a2, r2):
        return d + a1*np.exp(-t/r1) + a2*np.exp(-t/r2)

    def show(self, cell=None):
        """
        Plot results from run_iv()
        """
        if not HAVE_PG:
            raise Exception("Requires pyqtgraph")
        
        #
        # Generate figure with subplots
        #
        app = pg.mkQApp()
        win = pg.GraphicsWindow('%s  %s (%s)' % (cell.status['name'], cell.status['type'], cell.status['species']))
        self.win = win
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
        IVplot.plot(Icmd, self.peak_vm(), symbol='o', symbolBrush=(50, 150, 50, 255))
        IVplot.plot(Icmd, self.steady_vm(), symbol='s')


        # F/I relationship and raster plot
        spikes = self.spike_times()
        for i,times in enumerate(spikes):
            spikePlot.plot(x=times, y=[Icmd[i]]*len(times), pen=None, 
                           symbol='d', symbolBrush=colors[i])
        FIplot.plot(x=Icmd, y=[len(s) for s in spikes], symbol='o')
        
        
        # Print Rm, Vrest 
        (s, i, tau) = self.input_resistance_tau()
        print ("\nMembrane resistance (chord): {0:0.1f} MOhm  Taum1: {1:0.2f}  Taum2: {2:0.2f}".format(s, tau[0], tau[1]))
        ivals = np.array([Icmd.min(), Icmd.max()])
        vvals = s * ivals + i
        line = pg.QtGui.QGraphicsLineItem(ivals[0], vvals[0], ivals[1], vvals[1])
        line.setPen(pg.mkPen(255, 0, 0, 70))
        line.setZValue(-10)
        IVplot.addItem(line, ignoreBounds=True)
        
        print "Resting membrane potential: %0.1f mV\n" % self.rest_vm()
