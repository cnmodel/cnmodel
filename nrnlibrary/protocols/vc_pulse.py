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
