import os
import os.path
from neuron import h
import pylibrary.Utility as U
import pylibrary.PlotHelpers as PH
import numpy as np
import scipy
import scipy.integrate
import scipy.stats

from .protocol import Protocol


try:
    import pyqtgraph as pg
    HAVE_PG = True
except ImportError:
    HAVE_PG = False

from ..util.stim import make_pulse

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

    def run(self, vcrange, cell, dt=0.025):
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
        self.durs = [vstim.dur1, vstim.dur2, vstim.dur3]
        self.amps = [vstim.amp1, vstim.amp2, vstim.amp3]
        self.voltage_cmd = []
        tend = 900.0
        iv_nstepv = int(np.ceil((vmax - vmin) / vstep))
        iv_minv = vmin
        iv_maxv = vmax
        vstep = (iv_maxv - iv_minv) / iv_nstepv
        for i in range(iv_nstepv):
            self.voltage_cmd.append(float(i * vstep) + iv_minv)
        nreps = iv_nstepv
        h.dt = dt
        self.dt = h.dt
        for i in range(nreps):
            # Connect recording vectors
            self['v_soma'] = cell.soma(0.5)._ref_v
            self['i_inj'] = vstim._ref_i
            self['time'] = h._ref_t
            vstim.amp2 = self.voltage_cmd[i]
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

    def peak_im(self, window=0.4):
        """
        :param window: fraction of window to use for peak measurement, taken
        immediately following the beginning of the step
        Return steady-state membrane current for each trace.
        """
        Im = self.current_traces
        steps = len(Im)
        peakStop = (self.durs[0] + window*self.durs[1]) / self.dt
        peakStart = self.durs[0] / self.dt
        Vhold = self.amps[0] # np.mean([self.voltage_traces[i][:peakStart].mean() for i in range(steps)])
        Ipeak = []
        for i in range(steps):
            if self.voltage_cmd[i] > Vhold:
                Ipeak.append(Im[i][peakStart:peakStop].max())
            else:
                Ipeak.append(Im[i][peakStart:peakStop].min())
        return np.array(Ipeak)

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
        if cell is not None:
            win = pg.GraphicsWindow('%s  %s (%s)' % (cell.status['name'], cell.status['type'], cell.status['species']))
        else:
            win = pg.GraphisWindow('Voltage Clamp')
        self.win = win
        win.resize(1000, 800)
        Iplot = win.addPlot(labels={'left': 'Im (nA)', 'bottom': 'Time (ms)'})
        rightGrid = win.addLayout(rowspan=2)
        win.nextRow()
        Vplot = win.addPlot(labels={'left': 'V (mV)', 'bottom': 'Time (ms)'})

        IVplot = rightGrid.addPlot(labels={'left': 'Vm (mV)', 'bottom': 'Icmd (nA)'})
        IVplot.showGrid(x=True, y=True)
        rightGrid.nextRow()

        win.ci.layout.setRowStretchFactor(0, 10)
        win.ci.layout.setRowStretchFactor(1, 5)

        #
        # Plot simulation and analysis results
        #
        Vm = self.voltage_traces
        Iinj = self.current_traces
        Vcmd = self.voltage_cmd
        t = self.time_values
        steps = len(Vcmd)

        # plot I, V traces
        colors = [(i, steps*3./2.) for i in range(steps)]
        for i in range(steps):
            Vplot.plot(t, Vm[i], pen=colors[i])
            Iplot.plot(t, Iinj[i], pen=colors[i])

        # I/V relationships
        IVplot.plot(Vcmd, self.peak_im(), symbol='o', symbolBrush=(50, 150, 50, 255))
        IVplot.plot(Vcmd, self.steady_im(), symbol='s')

