from scipy import interpolate
import numpy as np
import matplotlib.pylab as mpl
import pyqtgraph as pg

from neuron import h

import nrnlibrary.util as util
from .protocol import Protocol
from .. import cells
from ..synapses import GluPSD, GlyPSD

mpl.rcParams['interactive'] = False



class SynapseTest(Protocol):
    def reset(self):
        super(SynapseTest, self).reset()

    def run(self, pre_cell, cell, n_synapses, temp=34.0):
        """ 
        Basic synapse test.
        Creates a presynaptic HH neuron and connects it to *cell* with 
        *n_synapses*.
        
        """
        synapses = []
        for i in range(n_synapses):
            synapses.append(pre_cell.connect(pre_cell.soma, cell.soma))
        
        self.synapses = synapses
        self.pre_cell = pre_cell
        self.allpsd = []
        # collect all PSDs across all synapses
        for syn in synapses:
            self.allpsd.extend(syn.psd.psd)

        VCLAMP = True
        glyPlot = False
        releasePlot = True
        
        #
        # voltage clamp the target cell
        #
        if VCLAMP == True:
            clampV = -65.0
            vccontrol = h.VClamp(0.5, sec=cell.soma)
            vccontrol.dur[0] = 10.0
            vccontrol.amp[0] = clampV
            vccontrol.dur[1] = 100.0
            vccontrol.amp[1] = clampV
            vccontrol.dur[2] = 20.0
            vccontrol.amp[2] = clampV

        #
        # set up stimulation of the presynaptic axon/terminal
        #
        istim = h.iStim(0.5, sec=pre_cell.soma)
        stim = {}
        stim['NP'] = 10
        stim['Sfreq'] = 100.0 # stimulus frequency
        stim['delay'] = 10.0
        stim['dur'] = 0.5
        stim['amp'] = 10.0
        stim['PT'] = 0.0
        stim['dt'] = h.dt
        (secmd, maxt, tstims) = util.make_pulse(stim)
        self.stim = stim
        
        istim.delay = 0
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0

        # istim current pulse train
        i_stim_vec = h.Vector(secmd)
        i_stim_vec.play(istim._ref_i, h.dt, 0)


        #
        # set up recordings
        #
        #self = {}
        #for var in ['v_pre', 'v_soma', 'i_soma', 'coh', 't', 'C0', 'C1', 'i_stim',
                    #'C2', 'C3', 'D', 'O1', 'O2', 'D1', 'D2', 'D3', 'Open', 'nmOpen', 'amOpen']:
            #self[var] = h.Vector()
        

        #for i in range(0, nANTerminals_ReleaseZones): 
            #self['isyn%03d' % i] = h.Vector(nANTerminals_ReleaseZones, 1000)
        
        # create hoc vectors for each parameter we wish to monitor and display
        synapse = synapses[0]
        self['v_pre'] = pre_cell.soma(0.5)._ref_v
        self['t'] = h._ref_t
        self['v_soma'] = pre_cell.soma(0.5)._ref_v
        self['coh'] = synapse.terminal.relsite._ref_XMTR[0]

        # make a synapse monitor for each release zone
        k = 0
        psd = self.allpsd
        for p in psd:
            #self['isyn%03d' % k] = h.Vector(len(psd), 1000)
            self['isyn%03d' % k] = psd[k]._ref_i
            k = k + 1
        
        self['Open'] = psd[0]._ref_Open
        if synapse.psd.kNMDA >= 0:
            self['nmOpen'] = psd[synapse.psd.kNMDA]._ref_Open
        if synapse.psd.kAMPA >= 0:
            self['amOpen'] = psd[synapse.psd.kAMPA]._ref_Open
        
        if isinstance(synapse.psd, GlyPSD):
            if synapse.psd.psdType == 'glyslow':
                nstate = 7
                self['C0'] = psd[0]._ref_C0
                self['C1'] = psd[0]._ref_C1
                self['C2'] = psd[0]._ref_C2
                self['O1'] = psd[0]._ref_O1
                self['O2'] = psd[0]._ref_O2
                self['D1'] = psd[0]._ref_D1
                #self['D3'] = psd[0]._ref_D3
                #self['O1'] = psd[0]._ref_O1
            elif synapse.psd.psdType == 'glyfast':
                nstate = 7
                self['C0'] = psd[0]._ref_C0
                self['C1'] = psd[0]._ref_C1
                self['C2'] = psd[0]._ref_C2
                self['C3'] = psd[0]._ref_C3
                self['O1'] = psd[0]._ref_O1
                self['O2'] = psd[0]._ref_O2

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
                isoma = np.zeros_like(self['isyn000'])
            for p in psd:
                thissyn = 'isyn%03d' % k
                if thissyn in self._vectors.keys():
                    isoma = isoma + self[thissyn]
                    k = k + 1
        self.isoma = isoma
        
        print 'Synapse.py: all runs done'
    

    def analyze(self, releasePlot=True, glyPlot=False):
        #
        # Analysis
        #
        synapse = self.synapses[0]
        nreq = 0
        nrel = 0
        nANTerminals = len(self.synapses)
        coh = [syn.terminal.relsite for syn in self.synapses]
        ntrel = np.zeros(nANTerminals)
        nANTerminals_ReleaseZones = synapse.terminal.n_rzones
        psd = self.allpsd
        
        #
        # compute some parameters
        #
        for j in range(0, nANTerminals):
            nreq = nreq + coh[j].nRequests # number of release requests during the for a terminal
            nrel = nrel + coh[j].nReleases # number of actual release events
            ntrel[j] = ntrel[j] + coh[j].nReleases # cumulative release events. (seems redundant)
            print 'Spikes: T%3d: = %3d ' % (j, coh[j].nRequests),
            print ' Releases = %4d from %d zones' % (coh[j].nReleases, nANTerminals_ReleaseZones)

        t = self['t']

        for i in range(0, nANTerminals):
            print 'ntrel[%d] = %d' % (i, ntrel[i])
        nreq = (nreq * nANTerminals_ReleaseZones)
        print 'Prel: %8.3f\n' % (coh[0].Dn * coh[0].Fn)
        print 'nreq: %d\n' % nreq
        if nreq > 0:
            print 'Rel Prob: %8.3f\n' % (float(nrel) / nreq)
        if synapse.psd.kNMDA >= 0:
            nmOmax = self['nmOpen'].max()
            amOmax = self['amOpen'].max()
            print 'Synapse.py: Max NMDAR Open Prob: %f   AMPA Open Prob: %f\n' % (nmOmax, amOmax)
            nmImax = self['isyn%03d' % synapse.psd.kNMDA].max()
            amImax = self['isyn%03d' % synapse.psd.kAMPA].max()
            if nmImax + amImax > 0.0:
                print 'Synapse.py: Max NMDAR I: %f   AMPA I: %f, N/(N+A): %f\n' % (
                    nmImax, amImax, nmImax / (nmImax + amImax))
            else:
                print "Synapse.py: release might have failed"

        #
        # plot the results for comparison
        #
        mpl.figure(1)
        g1 = mpl.subplot2grid((5, 1), (0, 0))
        p1 = g1.plot(t, self['v_pre'], color='black')
        g1.axes.set_ylabel('V pre')
        
        plotFocus = 'EPSC'
        
        if plotFocus == 'EPSC':
            g2 = mpl.subplot2grid((5, 1), (1, 0), rowspan=4)
            g2.plot(t, self.isoma, color='red')
            g2.axes.set_ylabel('I post')
            g2.axes.set_xlabel('Time (ms)')
        else:
            g2 = mpl.subplot2grid((5, 1), (1, 0), rowspan=1)
            g2.plot(t, self.isoma, color='cyan')
            g3 = mpl.subplot2grid((5, 1), (2, 0))
            g3.plot(t, self['v_pre'], color='blue')
            g3.plot(t, self['v_soma'], color='red')
            g4 = mpl.subplot2grid((5, 1), (3, 0))
            p4 = g4.plot(t, self['coh']) # glutamate
            g4.axes.set_ylabel('coh')
            g5 = mpl.subplot2grid((5, 1), (4, 0))
            k = 0
            for p in self.synapse.psd:
                if p.hname().find('NMDA', 0, 6) >= 0:
                    g5.plot(t, self['isyn%03d' % synapse.kNMDA]) # current through nmdar
                k = k + 1
            g5.axes.set_ylabel('inmda')
            g6 = mpl.subplot2grid((5, 1), (5, 0))
            k = 0
            for p in self.synapse.psd:
                if p.hname().find('NMDA', 0, 6) < 0:
                    g6.plot(t, self['isyn%03d' % synapse.kAMPA]) # glutamate
                k = k + 1
            g6.axes.set_ylabel('iAMPA')

        # Analyze the individual events. EPSCs get rise time, latency, half-width, and decay tau estimates.
        stim = self.stim
        ipi = 1000.0 / stim['Sfreq'] # convert from Hz (seconds) to msec.
        textend = 0.25 # allow response detection into the next frame
        pscpts = int((ipi + textend) / h.dt)
        tpsc = np.arange(0, ipi + textend, h.dt)
        ipsc = np.zeros((stim['NP'], pscpts))
        mpl.figure(num=220, facecolor='w')
        gpsc = mpl.subplot2grid((5, 2), (0, 0), rowspan=2, colspan=2)
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
            ipsc[i, :] = -self.isoma[istart:iend]
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
                    mpl.plot(t, self[var])
                    mpl.ylabel(var)
                    i = i + 1
            if psdType == 'glyfast':
                mpl.figure(2)
                for var in ['C0', 'C1', 'C2', 'C3', 'O1', 'O2', 'Open']:
                    mpl.subplot(7, 1, i + 1)
                    mpl.plot(t, self[var])
                    mpl.ylabel(var)
                    i = i + 1
        mpl.draw()
        mpl.show()
