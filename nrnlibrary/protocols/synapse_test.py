from scipy import interpolate
import numpy as np
import matplotlib.pylab as mpl
import pyqtgraph as pg

from neuron import h

from nrnlibrary.synapses import (Params, stochastic_synapses)
import nrnlibrary.util as util
from .protocol import Protocol
from .. import cells

mpl.rcParams['interactive'] = False



class SynapseTest(Protocol):
    def reset(self):
        super(IVCurve, self).reset()

    def run(self, cell, synapse, temp=34.0):
        """ 
        Test the synapse function. 
        Creates a presynaptic HH neuron and connects it to *cell* via *synapse*.
        
        v_pre is the presynaptic voltage
            v_soma is the postsynaptic voltage
            v_calyx is the calyx voltage
            coh is the calyx of held state
            isyn is the synaptic current
        """
        #
        # create presynaptic cell and wire up network
        #
        pre_cell = cells.HH()
        synapse.connect(pre_cell.soma, cell.soma)


        VCLAMP = True
        glyPlot = False
        releasePlot = True
        
        #
        # voltage clamp the target cell
        #
        if VCLAMP == True:
            clampV = -65.0
            vccontrol = h.VClamp(0.5, sec=SOMA)
            vccontrol.dur[0] = 10.0
            vccontrol.amp[0] = clampV
            vccontrol.dur[1] = 100.0
            vccontrol.amp[1] = clampV
            vccontrol.dur[2] = 20.0
            vccontrol.amp[2] = clampV

        #
        # adjust NMDA receptors
        #
        k = 0
        kNMDA = -1
        kAMPA = -1
        for p in psd:
            if p.hname().find('NMDA', 0, 6) >= 0:
                if TargetCellName == 'tstellate':
                    p.gNAR = 1.28 * AN_Po_Ratio * NMDARatio # for T-stellate cells, this yields correct AMPA, NMDA ratio of 1.13 at +40 mV
                if TargetCellName == 'bushy':
                    p.gNAR = 0.36 * AN_Po_Ratio * NMDARatio # for bushy cells, this yields correct AMPA, NMDA ratio of 0.429 at +40 mV
                    #if p is psd[0]:
                    #    print "NMDAR's for bushy cells set to %8.3f" % p.gNAR
                if TargetCellName == 'dstellate':
                    p.gNAR = 1.28 * AN_Po_Ratio * NMDARatio # same as for T-stellate (no other data)
                p.vshift = 0
                if kNMDA == -1:
                    kNMDA = k # save the first instance where we have an NMDA receptor
            else:
                if kAMPA == -1: # not NMDA, so get AMPA 
                    kAMPA = k
            k = k + 1
        
        #
        # set up stimulation of the presynaptic axon/terminal
        #
        istim = h.iStim(0.5, sec=axon[0])
        stim = {}
        stim['NP'] = 10
        stim['Sfreq'] = 100.0 # stimulus frequency
        stim['delay'] = 10.0
        stim['dur'] = 0.5
        stim['amp'] = 10.0
        stim['PT'] = 0.0
        stim['dt'] = h.dt
        (secmd, maxt, tstims) = util.make_pulse(stim)
        
        istim.delay = 0
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0

        #
        # set up recordings
        #
        vec = {}
        for var in ['v_pre', 'v_soma', 'i_soma', 'v_calyx', 'coh', 't', 'C0', 'C1', 'i_stim',
                    'C2', 'C3', 'D', 'O1', 'O2', 'D1', 'D2', 'D3', 'Open', 'nmOpen', 'amOpen']:
            vec[var] = h.Vector()
        
        # istim current pulse train
        vec['i_stim'] = h.Vector(secmd)

        # make a synapse monitor for each release zone
        for i in range(0, nANTerminals_ReleaseZones): 
            vec['isyn%03d' % i] = h.Vector(nANTerminals_ReleaseZones, 1000)
        
        # create hoc vectors for each parameter we wish to monitor and display
        vec['v_pre'].record(axon[0](0.5)._ref_v, sec=axon[0])
        vec['v_calyx'].record(calyx[0](0.5)._ref_v, sec=calyx[0])
        vec['t'].record(h._ref_t)
        vec['v_soma'].record(SOMA(0.5)._ref_v, sec=SOMA)
        vec['coh'].record(coh[0]._ref_XMTR[0], sec=TC)
        vec['i_stim'].play(istim._ref_i, h.dt, 0, sec=axon[0])

        k = 0
        for p in psd:
            vec['isyn%03d' % k] = h.Vector(len(psd), 1000)
            vec['isyn%03d' % k].record(psd[k]._ref_i, sec=TC)
            k = k + 1
        vec['Open'].record(psd[0]._ref_Open, sec=TC)
        if kNMDA >= 0:
            vec['nmOpen'].record(psd[kNMDA]._ref_Open, sec=TC)
        if kAMPA >= 0:
            vec['amOpen'].record(psd[kAMPA]._ref_Open, sec=TC)
        if psdType == 'glyslow':
            nstate = 7
            vec['C0'].record(psd[0]._ref_C0, sec=TC)
            vec['C1'].record(psd[0]._ref_C1, sec=TC)
            vec['C2'].record(psd[0]._ref_C2, sec=TC)
            vec['O1'].record(psd[0]._ref_O1, sec=TC)
            vec['O2'].record(psd[0]._ref_O2, sec=TC)
            vec['D1'].record(psd[0]._ref_D1, sec=TC)
            #vec['D3'].record(psd[0]._ref_D3, sec=TC)
            #vec['O1'].record(psd[0]._ref_O1, sec=TC)
        if psdType == 'glyfast':
            nstate = 7
            vec['C0'].record(psd[0]._ref_C0, sec=TC)
            vec['C1'].record(psd[0]._ref_C1, sec=TC)
            vec['C2'].record(psd[0]._ref_C2, sec=TC)
            vec['C3'].record(psd[0]._ref_C3, sec=TC)
            vec['O1'].record(psd[0]._ref_O1, sec=TC)
            vec['O2'].record(psd[0]._ref_O2, sec=TC)

        #
        # Run simulation
        #
        h.tstop = 200.0 # duration of a run
        h.celsius = temp
        tvec = np.arange(0, h.tstop, h.dt)
        h.init() # set up the parameters
        for nrep in xrange(1): # could do multiple runs.... 
            h.run()
            k = 0
            if nrep is 0: # save the soma current
                isoma = np.zeros_like(vec['isyn000'].to_python())
            for p in psd:
                thissyn = 'isyn%03d' % k
                if thissyn in vec.keys():
                    isoma = isoma + vec[thissyn].to_python()
                    k = k + 1
        print 'Synapse.py: all runs done'
    
    def analyze(self):
        #
        # Analysis
        #
        nreq = 0
        nrel = 0
        ntrel = np.zeros(nANTerminals)
        # compute some parameters
        for j in range(0, nANTerminals):
            nreq = nreq + coh[j].nRequests # number of release requests during the for a terminal
            nrel = nrel + coh[j].nReleases # number of actual release events
            ntrel[j] = ntrel[j] + coh[j].nReleases # cumulative release events. (seems redundant)
            print 'Spikes: T%3d: = %3d ' % (j, coh[j].nRequests),
            print ' Releases = %4d from %d zones' % (coh[j].nReleases, nANTerminals_ReleaseZones)

        for i in range(0, nANTerminals):
            print 'ntrel[%d] = %d' % (i, ntrel[i])
        nreq = (nreq * nANTerminals_ReleaseZones)
        print 'Prel: %8.3f\n' % (coh[0].Dn * coh[0].Fn)
        print 'nreq: %d\n' % nreq
        if nreq > 0:
            print 'Rel Prob: %8.3f\n' % (float(nrel) / nreq)
        if kNMDA >= 0:
            nmOmax = np.asarray(vec['nmOpen']).max()
            amOmax = np.asarray(vec['amOpen']).max()
            print 'Synapse.py: Max NMDAR Open Prob: %f   AMPA Open Prob: %f\n' % (nmOmax, amOmax)
            nmImax = np.asarray(vec['isyn%03d' % kNMDA]).max()
            amImax = np.asarray(vec['isyn%03d' % kAMPA]).max()
            if nmImax + amImax > 0.0:
                print 'Synapse.py: Max NMDAR I: %f   AMPA I: %f, N/(N+A): %f\n' % (
                    nmImax, amImax, nmImax / (nmImax + amImax))
            else:
                print "Synapse.py: release might have failed"

        # plot the results for comparison
        mpl.figure(1)
        g1 = mpl.subplot2grid((5, 1), (0, 0))
        p1 = g1.plot(vec['t'], vec['v_pre'], color='black')
        g1.axes.set_ylabel('V pre')
        if plotFocus == 'EPSC':
            g2 = mpl.subplot2grid((5, 1), (1, 0), rowspan=4)
            g2.plot(vec['t'].to_python(), isoma, color='red')
            g2.axes.set_ylabel('I post')
            g2.axes.set_xlabel('Time (ms)')
        else:
            g2 = mpl.subplot2grid((5, 1), (1, 0), rowspan=1)
            g2.plot(vec['t'].to_python(), isoma, color='cyan')
            g3 = mpl.subplot2grid((5, 1), (2, 0))
            g3.plot(vec['t'], vec['v_calyx'], color='blue')
            g3.plot(vec['t'], vec['v_soma'], color='red')
            g4 = mpl.subplot2grid((5, 1), (3, 0))
            p4 = g4.plot(vec['t'], vec['coh']) # glutamate
            g4.axes.set_ylabel('coh')
            g5 = mpl.subplot2grid((5, 1), (4, 0))
            k = 0
            for p in psd:
                if p.hname().find('NMDA', 0, 6) >= 0:
                    g5.plot(vec['t'], vec['isyn%03d' % kNMDA]) # current through nmdar
                k = k + 1
            g5.axes.set_ylabel('inmda')
            g6 = mpl.subplot2grid((5, 1), (5, 0))
            k = 0
            for p in psd:
                if p.hname().find('NMDA', 0, 6) < 0:
                    g6.plot(vec['t'], vec['isyn%03d' % kAMPA]) # glutamate
                k = k + 1
            g6.axes.set_ylabel('iAMPA')

        # Analyze the individual events. EPSCs get rise time, latency, half-width, and decay tau estimates.
        ipi = 1000.0 / stim['Sfreq'] # convert from Hz (seconds) to msec.
        textend = 0.25 # allow response detection into the next frame
        pscpts = int((ipi + textend) / h.dt)
        tpsc = np.arange(0, ipi + textend, h.dt)
        ipsc = np.zeros((stim['NP'], pscpts))
        mpl.figure(num=220, facecolor='w')
        gpsc = mpl.subplot2grid((5, 2), (0, 0), rowspan=2, colspan=2)
        ti = np.asarray(vec['t'])
        psc_20_lat = np.zeros((stim['NP'], 1)) # latency to 20% of rising amplitude
        psc_80_lat = np.zeros((stim['NP'], 1)) # latency to 80% of rising amplitude
        psc_hw = np.zeros((stim['NP'], 1)) # width at half-height
        psc_rt = np.zeros((stim['NP'], 1)) # 20-80 rise time
        tp = np.zeros((stim['NP'], 1))
        minLat = 0.5 # minimum latency for an event, in ms
        #    print 'NP: ', stim['NP']
        for i in range(stim['NP']):
            tstart = stim['delay'] + i * ipi
            istart = int(tstart / h.dt)
            minStart = int(minLat / h.dt)
            tp[i] = tstart - stim['delay']
            iend = int(istart + ((ipi + textend) / h.dt))
            #        print 'istart: %d iend: %d, len(isoma): %d\n' % (istart, iend, len(isoma))
            ipsc[i, :] = -isoma[istart:iend]
            psc_pk = np.argmax(ipsc[i, minStart:]) # position of the peak
            psc_pk = psc_pk + minStart - 1
            print 'i, pscpk, ipsc[i,pscpk]: ', i, psc_pk, ipsc[i, psc_pk]
            #       print 'minLat: %f   ipi+textend: %f, hdt: %f' % ((minLat, ipi+textend, h.dt))
            if psc_pk == 0:
                continue
            pkval = ipsc[i, psc_pk]
            psc_20_lat[i] = util.find_point(tpsc, ipsc[i, :], psc_pk, 0.2, direction='left', 
                                            limits=(minLat, ipi + textend, h.dt))
            psc_80_lat[i] = util.find_point(tpsc, ipsc[i, :], psc_pk, 0.8, direction='left', 
                                            limits=(minLat, ipi + textend, h.dt))
            psc_50l = util.find_point(tpsc, ipsc[i, :], psc_pk, 0.5, direction='left', 
                                    limits=(minLat, ipi + textend, h.dt))
            psc_50r = util.find_point(tpsc, ipsc[i, :], psc_pk, 0.5, direction='right', 
                                    limits=(minLat, ipi + textend, h.dt))
            if not np.isnan(psc_20_lat[i]) and not np.isnan(psc_80_lat[i]):
                psc_rt[i] = psc_80_lat[i] - psc_20_lat[i]
            else:
                psc_rt[i] = np.nan
            if not np.isnan(psc_50r) and not np.isnan(psc_50l):
                psc_hw[i] = float(psc_50r) - float(psc_50l)
                gpsc.plot(psc_50l, pkval * 0.5, 'k+')
                gpsc.plot(psc_50r, pkval * 0.5, 'k+')
                gpsc.plot(tpsc, ipsc[i, :].T)
            else:
                psc_hw[i] = np.nan
                gpsc.plot(tpsc, ipsc[i, :].T, color='0.6')
            gpsc.hold(True)
            gpsc.plot(psc_20_lat[i], pkval * 0.2, 'bo')
            gpsc.plot(psc_80_lat[i], pkval * 0.8, 'go')
            gpsc.plot(tpsc[psc_pk], pkval, 'ro')
        glat = mpl.subplot2grid((5, 2), (2, 0), colspan=2)
        grt = mpl.subplot2grid((5, 2), (3, 0), colspan=2)
        ghw = mpl.subplot2grid((5, 2), (4, 0), colspan=2)
        glat.plot(tp, psc_20_lat, 'bo-')
        glat.axes.set_ylabel('20%% Latency (ms)')
        grt.plot(tp, psc_rt, 'ms-')
        ghw.axes.set_ylabel('20-80 RT (ms)')
        ghw.plot(tp, psc_hw, 'rx-')
        ghw.axes.set_ylabel('Half-width (ms)')
        ghw.axes.set_xlabel('Time (ms)')

        mpl.show()
        #
        # now print some average values
        #
        nst = range(stim['NP'])
        analysisWindow = [nst[0:2], nst[-10:-1]]
        print analysisWindow
        print psc_rt
        RT_mean2080_early = np.ma.masked_invalid(psc_rt[analysisWindow[0]]).mean()
        RT_mean2080_late = np.ma.masked_invalid(psc_rt[analysisWindow[1]]).mean()
        Lat_mean20_early = np.ma.masked_invalid(psc_20_lat[analysisWindow[0]]).mean()
        Lat_mean20_late = np.ma.masked_invalid(psc_20_lat[analysisWindow[1]]).mean()
        HW_mean_early = np.ma.masked_invalid(psc_hw[analysisWindow[0]]).mean()
        HW_mean_late = np.ma.masked_invalid(psc_hw[analysisWindow[1]]).mean()
        print "Means: --------------"
        print RT_mean2080_early
        print Lat_mean20_early
        print HW_mean_early
        print 'Early:   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_mean2080_early, Lat_mean20_early,
                                                                                HW_mean_early)
        print 'Late :   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_mean2080_late, Lat_mean20_late,
                                                                                HW_mean_late)
        RT_std2080_early = np.ma.masked_invalid(psc_rt[analysisWindow[0]]).std()
        RT_std2080_late = np.ma.masked_invalid(psc_rt[analysisWindow[1]]).std()
        Lat_std20_early = np.ma.masked_invalid(psc_20_lat[analysisWindow[0]]).std()
        Lat_std20_late = np.ma.masked_invalid(psc_20_lat[analysisWindow[1]]).std()
        HW_std_early = np.ma.masked_invalid(psc_hw[analysisWindow[0]]).std()
        HW_std_late = np.ma.masked_invalid(psc_hw[analysisWindow[1]]).std()
        print "Standard Deviations: --------------"
        print 'Early:   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_std2080_early, Lat_std20_early,
                                                                                HW_std_early)
        print 'Late :   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_std2080_late, Lat_std20_late,
                                                                                HW_std_late)
        print "-----------------"

        if releasePlot: # plot the release event distributions over time in a separate window
            fig_ev = mpl.figure(num=301, facecolor='w')
            ph1 = fig_ev.add_subplot(211)
            ph2 = fig_ev.add_subplot(212)
            nbins = 50
            for j in range(0, nANTerminals):
                nev = coh[j].ev_index
                dis = np.array(coh[j].EventDist)[0:nev] # actual events, not original distribution
                tim = np.array(coh[j].EventTime)[0:nev]
                #            print "dis[%03d]: nev=%d " % (j, nev),
                #            print "# of spike Requests detected: %d \n" % (coh[j].nRequests)
                #print dis
                (hist, binedges) = np.histogram(dis, 25, range=(0.0, 5.0))
                ph1.hist(dis, nbins, range=(0.0, 5.0))
                ph2.plot(tim, dis, 'ro')
                ph2.axes.set_ylim((0, 5.0))
            ph2.axes.set_ylim((0, 5.0))
            ph2.axes.set_ylabel('# release events')
            ph2.axes.set_xlabel('Time post stimulus (ms)')
            ph1.axes.set_ylabel('Event latency (ms)')
            ph1.axes.set_xlabel('Time into train (ms)')
            mpl.draw()

        i = 0
        if glyPlot:
            if psdType == 'glyslow':
                mpl.figure(2)
                for var in ['C0', 'C1', 'C2', 'O1', 'O1', 'D1', 'Open']:
                    mpl.subplot(nstate, 1, i + 1)
                    mpl.plot(vec['t'], vec[var])
                    mpl.ylabel(var)
                    i = i + 1
            if psdType == 'glyfast':
                mpl.figure(2)
                for var in ['C0', 'C1', 'C2', 'C3', 'O1', 'O2', 'Open']:
                    mpl.subplot(7, 1, i + 1)
                    mpl.plot(vec['t'], vec[var])
                    mpl.ylabel(var)
                    i = i + 1
        mpl.draw()
        mpl.show()
