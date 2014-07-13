import os
import os.path
from neuron import h
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
