from neuron import h
import numpy as np
import scipy
import scipy.integrate
import scipy.stats
import scipy.optimize

try:
    import pyqtgraph as pg
    HAVE_PG = True
except ImportError:
    HAVE_PG = False

from ..util.stim import make_pulse
from ..util import fitting
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
        ivrange : list of tuples
            Each item in the list is (min, max, step) describing a range of 
            levels to test. Range values are inclusive, so the max value may
            appear in the test values. Using multiple ranges allows finer 
            measurements in some ranges. 
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

        # Calculate current pulse levels
        icur = []
        if isinstance(ivrange, tuple):
            ivrange = [ivrange]
        for c in ivrange:  # unpack current levels
            try:
                (imin, imax, istep) = c # unpack a tuple... or list
            except:
                raise TypeError("run_iv arguments must be a list of tuples [(imin, imax, istep), ...]")
            nstep = np.floor((imax-imin)/istep) + 1
            icur.extend(imin + istep * np.arange(nstep))  # build the list

        self.current_cmd = np.array(sorted(icur))
        nsteps = self.current_cmd.shape[0]
        # Configure IClamp
        if durs is None:
            durs = [10.0, 100.0, 50.0]
            
        self.durs = durs
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
        
        self.dt = dt
        self.temp = temp
        vec = {}

        for i in range(nsteps):
            # Generate current command for this level 
            stim['amp'] = self.current_cmd[i]
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
            cell.vm0 = None
            cell.cell_initialize()  # initialize the cell to it's rmp
            if i == 0:
                self.custom_init()
                r = cell.compute_rmrintau(auto_initialize=False)
                Rin, tau, v = r['Rin'], r['tau'], r['v']
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
    
    def input_resistance_tau(self, vmin=-10, imax=0, return_fits=False):
        """
        Estimate resting input resistance and time constant.
        
        Parameters
        ----------
        vmin : float
            minimum voltage to use in computation relative to resting
        imax : float
            maximum current to use in computation.
        return_eval : bool
            If True, return lmfit.ModelFit instances for the subthreshold trace
            fits as well.
            
        Returns
        -------
        dict :
            Dict containing:
            * 'slope' and 'intercept' keys giving linear 
              regression for subthreshold traces near rest
            * 'tau' giving the average first-exponential fit time constant
            * 'fits' giving a record array of exponential fit data to subthreshold
              traces.
        
        Analyzes only traces hyperpolarizing pulse traces near rest, with no 
        spikes.
        """
        Vss = self.steady_vm()
        vmin += self.rest_vm()
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
            print("WARNING: Not enough traces to do linear reggression in "
                  "IVCurve.input_resistance_tau().")
            print('{0:<15s}: {1:s}'.format('vss', ', '.join(['{:.2f}'.format(v) for v in Vss])))
            print('{0:<15s}: {1:s}'.format('vmask', repr(vmask.astype(int))))
            print('{0:<15s}: {1:s} '.format('imask', repr(imask.astype(int))))
            print('{0:<15s}: {1:s}'.format('spikemask', repr(smask.astype(int))))
            raise Exception("Not enough traces to do linear regression (see info above).")
        
        # Use these to measure input resistance by linear regression.
        reg = scipy.stats.linregress(Icmd[mask], Vss[mask])
        (slope, intercept, r, p, stderr) = reg

        # also measure the tau in the same traces:
        pulse_start = int(self.durs[0] / self.dt)
        pulse_stop = int((self.durs[0] + self.durs[1]) / self.dt)

        fits = []
        fit_inds = []
        tx = self.time_values[pulse_start:pulse_stop].copy()
        for i, m in enumerate(mask):
            if not m or (self.rest_vm() - Vss[i]) <= 1:
                continue
            
            trace = self.voltage_traces[i][pulse_start:pulse_stop]
            
            # find first minimum in the trace
            min_ind = np.argmin(trace)
            min_val = trace[min_ind]
            min_diff = trace[0] - min_val
            tau_est = min_ind * self.dt * (1 - 1 / np.e)
            
            # Fit cell charging to single exponential
            fit = fitting.Exp1().fit(trace[:min_ind],
                                     x=tx[:min_ind],
                                     xoffset=(tx[0], 'fixed'),
                                     yoffset=(min_val, -120, 0),
                                     amp=(min_diff, 0, 50),
                                     tau=(tau_est, 0.1, 50))

            # find first maximum in the trace (following with first minimum)
            max_ind = np.argmax(trace[min_ind:]) + min_ind
            max_val = trace[max_ind]
            max_diff = min_val - max_val
            tau2_est = max_ind * self.dt * (1 - 1 / np.e)
            amp1_est = fit.params['amp'].value
            tau1_est = fit.params['tau'].value
            amp2_est = fit.params['yoffset'] - max_val
            
            # fit up to first maximum with double exponential, using prior
            # fit as seed.
            fit = fitting.Exp2().fit(trace[:max_ind],
                                     method='nelder',
                                     x=tx[:max_ind],
                                     xoffset=(tx[0], 'fixed'),
                                     yoffset=(max_val, -120, 0),
                                     amp1=(amp1_est, 0, 200),
                                     tau1=(tau1_est, 0.1, 50),
                                     amp2=(amp2_est, -200, 0),
                                     tau_ratio=(tau2_est/tau1_est, 2, 50),
                                     tau2='tau_ratio * tau1'
                                     )
            
            fits.append(fit)
            fit_inds.append(i)

        # convert fits to record array
        dtype = [(k, float) for k in fits[0].params] + [('index', int)]
        fit_data = np.empty(len(fits), dtype=dtype)
        for i, fit in enumerate(fits):
            for k,v in fit.params.items():
                fit_data[i][k] = v.value
            fit_data[i]['index'] = fit_inds[i]
        
        if 'tau' in fit_data.dtype.fields:
            tau = fit_data['tau'].mean()
        else:
            tau = fit_data['tau1'].mean()
        
        ret = {'slope': slope, 
                'intercept': intercept,
                'tau': tau,
                'fits': fit_data}
        
        if return_fits:
            return ret, fits
        else:
            return ret

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
        IVplot.plot(Icmd, self.peak_vm(), symbol='o', symbolBrush=(50, 150, 50, 255), symbolSize=4.0)
        IVplot.plot(Icmd, self.steady_vm(), symbol='s', symbolSize=4.0)


        # F/I relationship and raster plot
        spikes = self.spike_times()
        for i,times in enumerate(spikes):
            spikePlot.plot(x=times, y=[Icmd[i]]*len(times), pen=None, 
                           symbol='d', symbolBrush=colors[i], symbolSize=4.0)
        FIplot.plot(x=Icmd, y=[len(s) for s in spikes], symbol='o', symbolSize=4.0)
        
        
        # Print Rm, Vrest 
        rmtau, fits = self.input_resistance_tau(return_fits=True)
        s = rmtau['slope']
        i = rmtau['intercept']
        #tau1 = rmtau['fits']['tau1'].mean()
        #tau2 = rmtau['fits']['tau2'].mean()
        #print ("\nMembrane resistance (chord): {0:0.1f} MOhm  Taum1: {1:0.2f}  Taum2: {2:0.2f}".format(s, tau1, tau2))
        
        # Plot linear subthreshold I/V relationship
        ivals = np.array([Icmd.min(), Icmd.max()])
        vvals = s * ivals + i
        line = pg.QtGui.QGraphicsLineItem(ivals[0], vvals[0], ivals[1], vvals[1])
        line.setPen(pg.mkPen(255, 0, 0, 70))
        line.setZValue(-10)
        IVplot.addItem(line, ignoreBounds=True)
        
        # plot exponential fits
        for fit in fits:
            t = np.linspace(self.durs[0], self.durs[0]+self.durs[1], 1000)
            y = fit.eval(x=t)
            Vplot.plot(t, y, pen={'color': (100, 100, 0), 'style': pg.QtCore.Qt.DashLine})

            # plot initial guess
            #y = fit.eval(x=t, **fit.init_params.valuesdict())
            #Vplot.plot(t, y, pen={'color': 'b', 'style': pg.QtCore.Qt.DashLine})
            
        print "Resting membrane potential: %0.1f mV\n" % self.rest_vm()
