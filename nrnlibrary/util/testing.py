import os
import os.path
from neuron import h
from neuron import *
from nrnlibrary.pynrnutilities import *
import pylibrary.Utility as U
import pylibrary.PlotHelpers as PH
import numpy
import scipy
import scipy.integrate
import scipy.stats as SStat

import matplotlib as MP # must call first... before pylag/pyplot or backends
MP.use('Qt4Agg')

import matplotlib.gridspec as GS
import mpl_toolkits.axes_grid1.inset_locator as INSETS
#import inset_axes, zoomed_inset_axes
import mpl_toolkits.axes_grid1.anchored_artists as ANCHOR
# import AnchoredSizeBar

stdFont = 'Arial'
#MP.use('pdf')
import  matplotlib.pyplot as pylab
pylab.rcParams['interactive'] = False
pylab.rcParams['mathtext.default'] = 'sf'
# next setting allows pdf font to be readable in Adobe Illustrator
pylab.rcParams['pdf.fonttype'] = 42
pylab.rcParams['figure.facecolor'] = 'white'

print MP.__version__


def make_pulse(stim, pulsetype="square"):
    """
    Generate a pulse train for current / voltage command. Returns a tuple:
    
    * w : stimulus waveform
    * maxt : duration of waveform
    * tstims : index of each pulse in the train
    
    Parameters:
    stim : dict
        Holds parameters that determine stimulus shape:
        
        * delay : time before first pulse
        * Sfreq : frequency of pulses 
        * dur : duration of one pulse
        * amp : pulse amplitude
        * PT : delay between end of train and test pulse (0 for no test)
        * NP : number of pulses
        * hold : holding level (optional)
    pulsetype : str
        'square' : square pulses  (default)
        'exp' : psg-like pulses
        
    """
    delay = int(numpy.floor(stim['delay'] / h.dt))
    ipi = int(numpy.floor((1000.0 / stim['Sfreq']) / h.dt))
    pdur = int(numpy.floor(stim['dur'] / h.dt))
    posttest = int(numpy.floor(stim['PT'] / h.dt))
    ndur = 5
    if stim['PT'] == 0:
        ndur = 1
    
    maxt = h.dt * (stim['delay'] + (ipi * (stim['NP'] + 3)) +
        posttest + pdur * ndur)
    
    hold = stim.get('hold', None)
    
    w = numpy.zeros(floor(maxt / h.dt))
    if hold is not None:
        w += hold
    
    #   make pulse
    tstims = [0] * int(stim['NP'])
    if pulsetype == 'square':
        for j in range(0, int(stim['NP'])):
            t = (delay + j * ipi) * h.dt
            w[delay + ipi * j:delay + (ipi * j) + pdur] = stim['amp']
            tstims[j] = delay + ipi * j
        if stim['PT'] > 0.0:
            send = delay + ipi * j
            for i in range(send + posttest, send + posttest + pdur):
                w[i] = stim['amp']

    if pulsetype == 'exp':
        for j in range(0, int(stim['NP'])):
            for i in range(0, len(w)):
                if delay + ipi * j + i < len(w):
                    w[delay + ipi * j + i] += (stim['amp'] *
                         (1.0 - exp(i / (pdur / -3.0))) *
                         exp(-1.0 * (i - (pdur / 3.0)) / pdur))
            tstims[j] = delay + ipi * j
        if stim['PT'] > 0.0:
            send = delay + ipi * j
            for i in range(send + posttest, len(w)):
                w[i] += (stim['amp'] *
                    (1.0 - exp(-1.0 * i / (pdur / 3.0))) *
                    exp(-1.0 * (i - (pdur / 3.0)) / pdur))
    
    return(w, maxt, tstims)


def run_iv(ivrange, cell, durs=None, sites=None, scales=None, reppulse=None):
    """
    Run a current-clamp I/V curve and display results.
    
    Parameters:
    ivrange : tuple
        (min, max, step)
    cell : Cell
        The Cell instance to test.
    durs : tuple
        durations of (pre, pulse, post) regions of the command
    sites : 
    scales : 
    reppulse : 
        stimulate with pulse train
    """
    try:
        (imin, imax, istep) = ivrange # unpack the tuple...
    except:
        raise TypeError("run_iv argument 1 must be a tuple (imin, imax, istep)")
    #print "min max step: ", imin, imax, istep
    
    if durs is None:
        durs = [10.0, 100.0, 50.0]
    icur = []
    if reppulse is None:
        istim = h.IClamp2(0.5, sec=cell.soma) # use our new iclamp method
        istim.dur[0] = durs[0]
        istim.amp[0] = 0
        istim.dur[1] = durs[1]
        istim.amp[1] = 0.0 #-70.00
        istim.dur[2] = durs[2]
        istim.amp[2] = 0.0 # 0.045
        istim.dur[3] = 0
        istim.amp[3] = 0
        istim.dur[4] = 0
        istim.amp[4] = 0
        tend = numpy.sum(durs)
    else:
        #
        # set up stimulation with a pulse train
        #
        istim = h.iStim(0.5, sec=cell.soma)
        stim = {}
        stim['NP'] = 10
        stim['Sfreq'] = 50.0 # stimulus frequency
        stim['delay'] = 10.0
        stim['dur'] = 2
        stim['amp'] = 1.0
        stim['PT'] = 0.0
        istim.delay = 0
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0
        (secmd, maxt, tstims) = make_pulse(stim)
        tend = maxt
        # istim current pulse train

    iv_nstepi = int(numpy.ceil((imax - imin) / istep))
    iv_mini = imin
    iv_maxi = imax
    nreps = iv_nstepi
    istep = (iv_maxi - iv_mini) / iv_nstepi
    iv_nstepi = iv_nstepi + 1
    for i in range(iv_nstepi):
        icur.append(float(i * istep) + iv_mini)
    nsteps = iv_nstepi
    vec = {}
    f1 = pylab.figure(1)
    p1 = pylab.subplot2grid((4, 1), (0, 0), rowspan=3)
    p2 = pylab.subplot2grid((4, 1), (3, 0), rowspan=1)
    #p3a = f1.add_subplot(6,1,6)
    #p3b = f1.add_subplot(6,1,5)
    f3 = pylab.figure(2)
    p3 = pylab.subplot2grid((2, 2), (0, 0), rowspan=1)
    p3.axes.set_ylabel(r'# spikes')
    p3.axes.set_xlabel(r'$I_{inj} (nA)$')
    p4 = pylab.subplot2grid((2, 2), (1, 0), rowspan=1)
    p4.axes.set_ylabel(r'Trial')
    p4.axes.set_xlabel(r'Time (ms)')
    p5 = pylab.subplot2grid((2, 2), (0, 1), rowspan=1)
    p5.axes.set_ylabel(r'V (mV)')
    p5.axes.set_xlabel(r'$I_{inj} (nA)$')
    p6 = pylab.subplot2grid((2, 2), (1, 1), rowspan=1)
    PH.cleanAxes([p1, p2, p3, p4, p5, p6])

    f4 = pylab.figure(3)
    p41 = pylab.subplot2grid((4, 1), (0, 0), rowspan=2)

#    if message is not None:
#        print 'meas: ', dir(measseg(0.5))

    clist = ['k-', 'r-', 'b-', 'y-', 'g-']
    slist = ['ko', 'rx', 'gx', 'bx', 'mx']
    splist = numpy.zeros(nsteps)
    meanVss = numpy.zeros(nsteps)
    meanIss = numpy.zeros(nsteps)
    minVpk = numpy.zeros(nsteps)
    for i in range(nsteps):
        for var in ['v_soma', 'i_inj', 'time', 'm', 'h', 'ah', 'bh', 'am',
                    'bm', 'gh', 'ik', 'ina', 'inat', 'i_stim']:
            vec[var] = h.Vector()
        if sites is not None:
            for j in range(len(sites)):
                vec['v_meas_%d' % (j)] = h.Vector()
        if not reppulse:
            istim.amp[1] = icur[i]
        else:
            stim['Amp'] = icur[i]
            (secmd, maxt, tstims) = make_pulse(stim)
            vec['i_stim'] = h.Vector(secmd)
#print "current: %f" % icur[i]
        h.tstop = tend
        vec['v_soma'].record(cell.soma(0.5)._ref_v, sec=cell.soma)
        vec['ik'].record(cell.soma(0.5)._ref_ik, sec=cell.soma)
        natFlag = False
        try:
            vec['inat'].record(cell.soma(0.5)._ref_inat, sec=cell.soma)
            natFlag = True
        except:
            vec['ina'].record(cell.soma(0.5)._ref_ina, sec=cell.soma)
            pass
        if sites is not None:
            for j in range(len(sites)):
                if sites[j] is not None:
                    print 'section %d : ' % (j),
                    print sites[j]
                    vec['v_meas_%d' % (j)].record(
                        sites[j](0.5)._ref_v, sec=sites[j])
        vec['i_inj'].record(istim._ref_i, sec=cell.soma)
        vec['gh'].record(cell.soma().ihvcn._ref_i, sec=cell.soma)
        vec['time'].record(h._ref_t)
        if reppulse is not None:
            vec['i_stim'].play(istim._ref_i, h.dt, 0, sec=cell.soma)

        # h.t = -200.
        # dtsav = h.dt
        # h.dt = 1e9
        # while h.t < 0:
        #     h.fadvance()
        # h.dt = dtsav
        # h.t = 0
        # h.fcurrent()
        # h.cvode.re_init()
        h.init()
        h.run()
        #tvec = arange(0, h.tstop, h.dt)
        p1.plot(vec['time'], vec['v_soma'], 'k') # soma is plotted in black...
        p1.axes.set_ylabel('V (mV)')
        ik = numpy.asarray(vec['ik'])
        ina = numpy.asarray(vec['ina'])
        if natFlag:
            if len(ina) == 0:
                ina = numpy.asarray(vec['inat'])
            else:
                ina = ina + numpy.asarray(vec['inat'])
        t = numpy.asarray(vec['time'])
        iQ = scipy.integrate.trapz(ik, t) # total charge at end of run
        iQKt = scipy.integrate.cumtrapz(ik, t, initial=0.0)
        # cumulative with trapezoidal integration
        iQNat = scipy.integrate.cumtrapz(ina, t, initial=0.0)
        p41.plot(t, iQKt, 'g')
        p41.plot(t, iQNat, 'r')
        PH.cleanAxes(p41)
#        PH.cleanAxes(p1)
        mwine = durs[0] + durs[1]
        mwins = mwine - 0.2 * durs[1]
        vsoma = numpy.asarray(vec['v_soma'])
        (meanVss[i], r2) = U.measure('mean', vec['time'], vsoma, mwins, mwine)
        (meanIss[i], r2) = U.measure('mean', vec['time'], vec['i_inj'],
                            mwins, mwine)
        (minVpk[i], r2) = U.measure('min', vec['time'], vsoma, durs[0],
                            durs[0] + 0.5 * durs[1])
        if sites is not None:
            for j in range(len(sites)):
                if sites[j] is not None:
                    p1.plot(vec['time'], numpy.asarray(
                            vec['v_meas_%d' % (j)]), clist[j])
        p2.plot(vec['time'], vec['i_inj'], 'k')
#        PH.cleanAxes(p2)
        p2.axes.set_ylabel(r'$I_{inj} (nA)$')
        p2.axes.set_xlabel(r'Time (ms)')
        #p3b.plot(vec['time'], vec['gh'])
        spli = findspikes(vec['time'], vec['v_soma'], -30.0)
        nsoma = i * numpy.ones(len(spli))
        splist[i] = len(spli)
        p4.plot(spli, nsoma, 'bo-')
        if sites is not None:
            for j in range(len(sites)):
                if sites[j] is not None:
                    splim = U.findspikes(vec['time'], numpy.asarray(
                            vec['v_meas_%d' % (j)]), -30.0)
                    nseg = i * numpy.ones(len(splim))
                    if len(splim) > 0 and len(nseg) > 0:
                        p2.plot(splim, nseg, slist[j])

        pylab.draw()

    ok1 = numpy.where(meanIss <= 0.0)[0].tolist()
    ok2 = numpy.where(meanVss >= -70.0)[0].tolist()
    ok3 = numpy.where(splist == 0)[0].tolist()
    ok = list(set(ok1).intersection(set(ok2)))
    #Linear regression using stats.linregress
    if len(ok) > 1: # need 2 points to make that line
        (a_s, b_s, r, tt, stderr)=SStat.linregress(meanIss[ok], meanVss[ok])
        print('Linear regression using stats.linregress')
        print('regression: slope=%.2f intercept=%.2f, std error= %.3f'
         % (a_s, b_s, stderr))
        print '  r: %.3f   p: %.3f' % (r, tt)

    p2.set_xlim(0, 160)
    p1.set_xlim(0, 160)
    if scales is not None:
        p3.set_xlim(scales[0], scales[2])
        p3.set_ylim(scales[4], scales[5])
        PH.crossAxes(p5, limits=scales[0:4], xyzero=scales[9])
        if scales[6] == 'offset':
            PH.nice_plot(p3, direction='outward')
    p5.plot(meanIss[ok3], meanVss[ok3], 'ko-')
    p5.plot(meanIss[ok1], minVpk[ok1], 'ks-')
    p3.plot(icur, splist, 'ro-')

    print 'I,Vss,Vpk,SpikesperSec'
    for i in range(nsteps):
        print '%8.4f,%8.3f,%8.3f,%8.2f' % (icur[i],
            meanVss[i], minVpk[i], splist[i])


    pylab.show()


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
    iv_nstepv = int(numpy.ceil((vmax - vmin) / vstep))
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
    meanVss = numpy.zeros(nreps)
    meanIss = numpy.zeros(nreps)
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
            vc = numpy.asarray(vec['Vsoma'])
            tc = numpy.asarray(vec['time'])
            
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
