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
