from scipy import interpolate
import numpy as np
import matplotlib.pylab as mpl
import pyqtgraph as pg

from neuron import h

import nrnlibrary.util as util
from .protocol import Protocol
from .. import cells
from ..synapses import GluPSD, GlyPSD

#mpl.rcParams['interactive'] = False


class SynapseTest(Protocol):
    def reset(self):
        super(SynapseTest, self).reset()

    def run(self, pre_sec, post_sec, n_synapses, temp=34.0, dt=0.025):
        """ 
        Basic synapse test. Connects sections of two cells with *n_synapses*.
        The cells are allowed to negotiate the details of the connecting 
        synapse. The presynaptic soma is then driven with a pulse train
        followed by a recovery pulse of varying delay.
        
        Analyses:
        
        * Distribution of PSG amplitude, kinetics, and latency
        * Synaptic depression / facilitation and recovery timecourses
        """
        pre_cell = cells.cell_from_section(pre_sec)
        post_cell = cells.cell_from_section(post_sec)
        synapses = []
        for i in range(n_synapses):
            synapses.append(pre_cell.connect(pre_sec, post_sec))
        
        self.synapses = synapses
        self.pre_sec = pre_sec
        self.post_sec = post_sec
        self.pre_cell = pre_cell
        self.post_cell = post_cell
        
        #
        # voltage clamp the target cell
        #
        clampV = 40.0
        vccontrol = h.VClamp(0.5, sec=post_cell.soma)
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
        stim['dt'] = dt
        (secmd, maxt, tstims) = util.make_pulse(stim)
        self.stim = stim
        
        istim.delay = 0
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0

        # istim current pulse train
        i_stim_vec = h.Vector(secmd)
        i_stim_vec.play(istim._ref_i, dt, 0)


        # create hoc vectors for each parameter we wish to monitor and display
        synapse = synapses[0]
        self['v_pre'] = pre_cell.soma(0.5)._ref_v
        self['t'] = h._ref_t
        self['v_soma'] = pre_cell.soma(0.5)._ref_v
        self['relsite_xmtr'] = synapse.terminal.relsite._ref_XMTR[0]


        # make a synapse monitor for each release zone
        #self.all_psd = []
        self.all_nmda = []
        self.all_ampa = []
        for syn in synapses:
            # collect all PSDs across all synapses
            #self.all_psd.extend(syn.psd.all_psd)
            self.all_nmda.extend(syn.psd.nmda_psd)
            self.all_ampa.extend(syn.psd.ampa_psd)
        
        #  Record current through all PSDs individually
        for k,p in enumerate(self.all_nmda):
            self['iNMDA%03d' % k] = p._ref_i
            self['opNMDA%03d' % k] = p._ref_Open
        for k,p in enumerate(self.all_ampa):
            self['iAMPA%03d' % k] = p._ref_i
            self['opAMPA%03d' % k] = p._ref_Open
        
        # Record sample AMPA and NMDA currents
        #self['iNMDA'] = synapse.psd.nmda_psd[0]._ref_i
        #self['iAMPA'] = synapse.psd.ampa_psd[0]._ref_i
        
        #self['Open'] = self.all_psd[0]._ref_Open
        #if isinstance(synapse.psd, GluPSD):
            #if len(synapse.psd.nmda_psd) > 0:
               #self['nmOpen'] = synapse.psd.nmda_psd[0]._ref_Open
            #if len(synapse.psd.ampa_psd) > 0:
               #self['amOpen'] = synapse.psd.ampa_psd[0]._ref_Open
        
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
        h.dt = dt
        self.temp = temp
        self.dt = dt
        
        self.custom_init()
        for nrep in xrange(1): # could do multiple runs.... 
            h.run()
            
            if nrep is 0: # add up psd current across all runs
                isoma = np.zeros_like(self['iAMPA000'])
            for k in range(len(self.all_ampa)):
                isoma += self['iAMPA%03d'%k]
            for k in range(len(self.all_nmda)):
                isoma += self['iNMDA%03d'%k]
        self.isoma = isoma

    def release_events(self):
        """
        Analyze results and return a dict of values related to terminal release 
        probability:
            
            n_zones: Array containing the number of release zones for each
                     synapse.
            n_requests: Array containing number of release requests for each 
                        synapse. Note for multi-zone synapses, a single 
                        presynaptic spike results in one release request _per_
                        zone.
            n_releases: Array containing actual number of releases for each 
                        synapse.
            tot_requests: The total number of release requests across all
                          release zones. 
            tot_releases: The total number of releases.
            release_p: Release probability computed as 
                       tot_releases / tot_requests
        """
        synapse = self.synapses[0]
        
        ret = {}
        #
        # Count spikes and releases for each terminal
        #
        ret['n_zones'] = np.array([syn.terminal.n_rzones for syn in self.synapses])
        ret['n_spikes'] = np.array([syn.terminal.relsite.nRequests for syn in self.synapses])
        ret['n_requests'] = ret['n_spikes'] * ret['n_zones']
        ret['n_releases'] = np.array([syn.terminal.relsite.nReleases for syn in self.synapses])

        #
        # Compute release probability
        #
        # total number of release requests
        ret['tot_requests'] = ret['n_requests'].sum()
        # total number of actual release events        
        ret['tot_releases'] = ret['n_releases'].sum() 
        
        if ret['tot_requests'] > 0:
            ret['release_p'] = float(ret['tot_releases']) / ret['tot_requests']
        else:
            ret['release_p'] = np.nan
        
        return ret

    def open_probability(self):
        """ 
        Analyze results and return a dict of values related to psd open 
        probability:
            
            op_nmda:
            op_ampa:
            imax_nmda:
            imax_ampa:
        """
        synapse = self.synapses[0]
        if isinstance(synapse.psd, GluPSD) and len(synapse.psd.nmda_psd) > 0:
            # find a psd with ampa and nmda currents
            nmImax = 0
            amImax = 0
            nmOmax = 0
            amOmax = 0
            #self.win.nextRow()
            for i in range(len(self.all_ampa)):
                nm = np.abs(self['iNMDA%03d'%i]).max()
                am = np.abs(self['iAMPA%03d'%i]).max()
                opnm = np.abs(self['opNMDA%03d'%i]).max()
                opam = np.abs(self['opAMPA%03d'%i]).max()
                #plt = self.win.addPlot()
                #plt.plot(self['iNMDA%03d'%i])
                #plt = self.win.addPlot()
                #plt.plot(self['iAMPA%03d'%i])
                #plt = self.win.addPlot()
                #plt.plot(self['opNMDA%03d'%i])
                #plt = self.win.addPlot()
                #plt.plot(self['opAMPA%03d'%i])
                #self.win.nextRow()
                if nm != 0 or am != 0:
                    nmImax = nm
                    amImax = am
                    nmOmax = opnm
                    amOmax = opam
                    break
            
            return {
                'op_nmda': nmOmax,
                'op_ampa': amOmax,
                'imax_nmda': nmImax,
                'imax_ampa': amImax,
            }

    def analyze_events(self):
        stim = self.stim
        ipi = 1000.0 / stim['Sfreq'] # convert from Hz (seconds) to msec.
        t_extend = 0.25 # allow response detection into the next frame
        extend_pts = int(t_extend / self.dt)
        
        pscpts = int(ipi / self.dt) + extend_pts # number of samples to analyze for each psc
        ipsc = np.zeros((stim['NP'], pscpts))  # storage for psc currents
        tpsc = np.arange(0, ipi + t_extend, self.dt) # time values corresponding to *ipsc*
        
        #mpl.figure(num=220, facecolor='w')
        #gpsc = mpl.subplot2grid((5, 2), (0, 0), rowspan=2, colspan=2)
        psc_20_lat = np.zeros((stim['NP'], 1)) # latency to 20% of rising amplitude
        psc_80_lat = np.zeros((stim['NP'], 1)) # latency to 80% of rising amplitude
        psc_hw = np.zeros((stim['NP'], 1)) # width at half-height
        psc_rt = np.zeros((stim['NP'], 1)) # 20-80 rise time
        tp = np.zeros((stim['NP'], 1))  # pulse time relative to first pulse
        events = np.zeros(stim['NP'], dtype=[
            ('20% latency', float),
            ('80% latency', float),
            ('half width', float),
            ('half left', float),
            ('half right', float),
            ('rise time', float),
            ('pulse time', float),
            ('peak', float),
            ('peak index', int),
        ])
        
        minLat = 0.0 # minimum latency for an event, in ms
        minStart = int(minLat / self.dt)  # first index relative to pulse to search for psc peak
        
        for i in range(stim['NP']):
            tstart = stim['delay'] + i * ipi # pulse start time 
            istart = int(tstart / self.dt)   # pulse start index
            tp[i] = tstart - stim['delay']
            iend = istart + pscpts
            #        print 'istart: %d iend: %d, len(isoma): %d\n' % (istart, iend, len(isoma))
            ipsc[i, :] = -self.isoma[istart:iend]
            psc_pk = minStart + np.argmax(ipsc[i, minStart:-extend_pts]) # position of the peak
            #print 'i, pscpk, ipsc[i,pscpk]: ', i, psc_pk, ipsc[i, psc_pk]
            #       print 'minLat: %f   ipi+t_extend: %f, hdt: %f' % ((minLat, ipi+t_extend, self.dt))
            if psc_pk == minStart:
                continue
            pkval = ipsc[i, psc_pk]
            events['peak'][i] = pkval
            events['peak index'][i] = psc_pk
            
            lat20 = util.find_point(tpsc, ipsc[i, :], psc_pk, 0.2, direction='left', 
                                            limits=(minLat, ipi + t_extend, self.dt))
            lat80 = util.find_point(tpsc, ipsc[i, :], psc_pk, 0.8, direction='left', 
                                            limits=(minLat, ipi + t_extend, self.dt))
            events['20% latency'][i] = lat20
            events['80% latency'][i] = lat80
            
            psc_50l = util.find_point(tpsc, ipsc[i, :], psc_pk, 0.5, direction='left', 
                                    limits=(minLat, ipi + t_extend, self.dt))
            psc_50r = util.find_point(tpsc, ipsc[i, :], psc_pk, 0.5, direction='right', 
                                    limits=(minLat, ipi + t_extend, self.dt))
            events['half left'] = psc_50l
            events['half right'] = psc_50r
            
            if not np.isnan(lat20) and not np.isnan(lat80):
                events['rise time'][i] = lat80 - lat20
            else:
                events['rise time'][i] = np.nan
            if not np.isnan(psc_50r) and not np.isnan(psc_50l):
                events['half width'][i] = float(psc_50r) - float(psc_50l)
                #gpsc.plot(psc_50l, pkval * 0.5, 'k+')
                #gpsc.plot(psc_50r, pkval * 0.5, 'k+')
                #gpsc.plot(tpsc, ipsc[i, :].T)
            else:
                events['half width'][i] = np.nan

        print events
        return events

    def show(self, releasePlot=True, glyPlot=False, plotFocus='EPSC'):
        self.win = pg.GraphicsWindow()
        synapse = self.synapses[0]
        
        #
        # Print parameters related to release probability
        #
        events = self.release_events()
        ns = len(self.synapses)
        for i in range(ns):
            v = (i, events['n_spikes'][i], events['n_zones'][i], events['n_releases'][i])
            print 'Synapse %d:  spikes: %d  zones: %d  releases: %d' % v
        print ""
        print 'Total release requests: %d' % events['tot_requests']
        print 'Total release events:   %d' % events['tot_releases']
        print 'Release probability: %8.3f' % events['release_p']
        prel_final = synapse.terminal.relsite.Dn * synapse.terminal.relsite.Fn
        print 'Final release probability (Dn * Fn): %8.3f' % prel_final


        #
        # Compute NMDA / AMPA open probability
        #
        print ""
        oprob = self.open_probability()
        print 'Max NMDAR Open Prob: %f   AMPA Open Prob: %f' % (oprob['op_nmda'], oprob['op_ampa'])
        nmImax, amImax = oprob['imax_nmda'], oprob['imax_ampa']
        print 'Max NMDAR I: %f   AMPA I: %f' % (nmImax, amImax)
        if nmImax + amImax != 0.0:
            print '   N/(N+A): %f\n' % (nmImax / (nmImax + amImax))
        else:
            print "   (no NMDA/AMPA current; release might have failed)"


        #
        # plot the results for comparison
        #
        t = self['t']

        #mpl.figure(1)
        #g1 = mpl.subplot2grid((6, 1), (0, 0))
        #p1 = g1.plot(t, self['v_pre'], color='black')
        #g1.axes.set_ylabel('V pre')
        #g1.set_title(self.pre_cell.status['name'])
        
        p1 = self.win.addPlot(title=self.pre_cell.status['name'])
        p1.setLabels(left='V pre (mV)', bottom='Time (ms)')
        p1.plot(t, self['v_pre'])
        
        if plotFocus == 'EPSC':
            #g2 = mpl.subplot2grid((6, 1), (1, 0), rowspan=4)
            #g2.plot(t, self.isoma, color='red')
            #g2.axes.set_ylabel('I post')
            #g2.axes.set_xlabel('Time (ms)')
            #g2.set_title(self.post_cell.status['name'])
            self.win.nextRow()
            p2 = self.win.addPlot(title=self.post_cell.status['name'])
            p2.plot(t, self.isoma, pen='r')
            p2.setLabels(left='Total PSD current (nA)', bottom='Time (ms)')
        else:
            # todo: resurrect this?
            g2 = mpl.subplot2grid((6, 1), (1, 0), rowspan=1)
            g2.plot(t, self.isoma, color='cyan')
            g3 = mpl.subplot2grid((6, 1), (2, 0))
            g3.plot(t, self['v_pre'], color='blue')
            g3.plot(t, self['v_soma'], color='red')
            g4 = mpl.subplot2grid((6, 1), (3, 0))
            p4 = g4.plot(t, self['relsite_xmtr']) # glutamate
            g4.axes.set_ylabel('relsite_xmtr')
            g5 = mpl.subplot2grid((6, 1), (4, 0))
            for k,p in enumerate(synapse.psd.all_psd):
                if p.hname().find('NMDA', 0, 6) >= 0:
                    g5.plot(t, self['isyn%03d' % k]) # current through nmdar
            g5.axes.set_ylabel('inmda')
            g6 = mpl.subplot2grid((6, 1), (5, 0))
            for k,p in enumerate(synapse.psd.all_psd):
                if p.hname().find('NMDA', 0, 6) < 0:
                    g6.plot(t, self['isyn%03d' % k]) # glutamate
            g6.axes.set_ylabel('iAMPA')


        # 
        # Analyze the individual events. 
        # EPSCs get rise time, latency, half-width, and decay tau estimates.
        #
        events = self.analyze_events()
        
        #for i in range(len(events)):
            #if not np.isnan(events['half width'][i]):
                #gpsc.plot(psc_50l, pkval * 0.5, 'k+')
                #gpsc.plot(psc_50r, pkval * 0.5, 'k+')
                #gpsc.plot(tpsc, ipsc[i, :].T)
            #else:
                #gpsc.plot(tpsc, ipsc[i, :].T, color='0.6')
            #gpsc.hold(True)
            #gpsc.plot(events['20% latency'][i], pkval * 0.2, 'bo')
            #gpsc.plot(events['80% latency'][i], pkval * 0.8, 'go')
            #gpsc.plot(tpsc[psc_pk], pkval, 'ro')
        self.win.nextRow()
        plt = self.win.addPlot(labels={'left': '20% Latency (ms)', 'bottom': 'Pulse Time (ms)'})
        plt.plot(events['20% latency'], pen=None, symbol='o')
        self.win.nextRow()
        plt = self.win.addPlot(labels={'left': '80% Latency (ms)', 'bottom': 'Pulse Time (ms)'})
        plt.plot(events['80% latency'], pen=None, symbol='o')
        #glat = mpl.subplot2grid((5, 2), (2, 0), colspan=2)
        #grt = mpl.subplot2grid((5, 2), (3, 0), colspan=2)
        #ghw = mpl.subplot2grid((5, 2), (4, 0), colspan=2)
        #glat.plot(tp, events['20% latency'], 'bo-')
        #glat.axes.set_ylabel('20%% Latency (ms)')
        #grt.plot(tp, events['rise time'], 'ms-')
        #ghw.axes.set_ylabel('20-80 RT (ms)')
        #ghw.plot(tp, events['half width'], 'rx-')
        #ghw.axes.set_ylabel('Half-width (ms)')
        #ghw.axes.set_xlabel('Time (ms)')

        #mpl.show()
        return
        
        #
        # now print some average values
        #
        nst = range(stim['NP'])
        analysisWindow = [nst[0:2], nst[-10:-1]]
        print analysisWindow
        print events['rise time']
        RT_mean2080_early = np.ma.masked_invalid(events['rise time'][analysisWindow[0]]).mean()
        RT_mean2080_late = np.ma.masked_invalid(events['rise time'][analysisWindow[1]]).mean()
        Lat_mean20_early = np.ma.masked_invalid(events['20% latency'][analysisWindow[0]]).mean()
        Lat_mean20_late = np.ma.masked_invalid(events['20% latency'][analysisWindow[1]]).mean()
        HW_mean_early = np.ma.masked_invalid(events['half width'][analysisWindow[0]]).mean()
        HW_mean_late = np.ma.masked_invalid(events['half width'][analysisWindow[1]]).mean()
        print "Means: --------------"
        print RT_mean2080_early
        print Lat_mean20_early
        print HW_mean_early
        print 'Early:   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_mean2080_early, Lat_mean20_early,
                                                                                HW_mean_early)
        print 'Late :   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_mean2080_late, Lat_mean20_late,
                                                                                HW_mean_late)
        RT_std2080_early = np.ma.masked_invalid(events['rise time'][analysisWindow[0]]).std()
        RT_std2080_late = np.ma.masked_invalid(events['rise time'][analysisWindow[1]]).std()
        Lat_std20_early = np.ma.masked_invalid(events['20% latency'][analysisWindow[0]]).std()
        Lat_std20_late = np.ma.masked_invalid(events['20% latency'][analysisWindow[1]]).std()
        HW_std_early = np.ma.masked_invalid(events['half width'][analysisWindow[0]]).std()
        HW_std_late = np.ma.masked_invalid(events['half width'][analysisWindow[1]]).std()
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
            for j in range(0, len(self.synapses)):
                relsite = self.synapses[j].terminal.relsite
                nev = relsite.ev_index
                dis = np.array(relsite.EventDist)[0:nev] # actual events, not original distribution
                tim = np.array(relsite.EventTime)[0:nev]
                #            print "dis[%03d]: nev=%d " % (j, nev),
                #            print "# of spike Requests detected: %d \n" % (relsites[j].nRequests)
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
            if synapse.psd.psdType == 'glyslow':
                mpl.figure(2)
                for var in ['C0', 'C1', 'C2', 'O1', 'O1', 'D1', 'Open']:
                    mpl.subplot(nstate, 1, i + 1)
                    mpl.plot(t, self[var])
                    mpl.ylabel(var)
                    i = i + 1
            if synapse.psd.psdType == 'glyfast':
                mpl.figure(2)
                for var in ['C0', 'C1', 'C2', 'C3', 'O1', 'O2', 'Open']:
                    mpl.subplot(7, 1, i + 1)
                    mpl.plot(t, self[var])
                    mpl.ylabel(var)
                    i = i + 1
        mpl.draw()
        mpl.show()
