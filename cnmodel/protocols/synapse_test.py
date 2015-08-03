from collections import OrderedDict
from scipy import interpolate
import numpy as np
import pyqtgraph as pg

from neuron import h

import cnmodel.util as util
from .protocol import Protocol
from .. import cells
from ..synapses import GluPSD, GlyPSD
from ..util.find_point import find_crossing


class SynapseTest(Protocol):
    def reset(self):
        super(SynapseTest, self).reset()

    def run(self, pre_sec, post_sec, n_synapses, temp=34.0, dt=0.025, 
            vclamp=40.0, iterations=1, tstop=200.0, stim_params=None, **kwds):
        """ 
        Basic synapse test. Connects sections of two cells with *n_synapses*.
        The cells are allowed to negotiate the details of the connecting 
        synapse. The presynaptic soma is then driven with a pulse train
        followed by a recovery pulse of varying delay.
        
        *stim_params* is an optional dictionary with keys 'NP', 'Sfreq', 'delay',
        'dur', 'amp'.
        
        Analyses:
        
        * Distribution of PSG amplitude, kinetics, and latency
        * Synaptic depression / facilitation and recovery timecourses
        """
        Protocol.run(self, **kwds)
        
        pre_cell = cells.cell_from_section(pre_sec)
        post_cell = cells.cell_from_section(post_sec)
        synapses = []
        for i in range(n_synapses):
            synapses.append(pre_cell.connect(post_cell))
        
        self.synapses = synapses
        self.pre_sec = synapses[0].terminal.section
        self.post_sec = synapses[0].psd.section
        self.pre_cell = pre_cell
        self.post_cell = post_cell
        
        #
        # voltage clamp the target cell
        #
        clampV = vclamp
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
        stim = {
            'NP': 10,
            'Sfreq': 100.0,
            'delay': 10.0,
            'dur': 0.5,
            'amp': 10.0,
            'PT': 0.0,
            'dt': dt,
        }
        if stim_params is not None:
            stim.update(stim_params)
        (secmd, maxt, tstims) = util.make_pulse(stim)
        self.stim = stim

        if tstop is None:
            tstop = len(secmd) * dt
        
        istim.delay = 0
        istim.dur = 1e9 # these actually do not matter...
        istim.iMax = 0.0

        # istim current pulse train
        i_stim_vec = h.Vector(secmd)
        i_stim_vec.play(istim._ref_i, dt, 0)

        # create hoc vectors for each parameter we wish to monitor and display
        synapse = synapses[0]

        self.all_psd = []
        for syn in synapses:
            # collect all PSDs across all synapses
            self.all_psd.extend(syn.psd.all_psd)
            
        #for i, cleft in enumerate(synapse.psd.clefts):
            #self['cleft_xmtr%d' % i] = cleft._ref_CXmtr
            #self['cleft_pre%d' % i] = cleft._ref_pre
            #self['cleft_xv%d' % i] = cleft._ref_XV
            #self['cleft_xc%d' % i] = cleft._ref_XC
            #self['cleft_xu%d' % i] = cleft._ref_XU

        #
        # Run simulation
        #
        h.tstop = tstop # duration of a run
        h.celsius = temp
        h.dt = dt
        self.temp = temp
        self.dt = dt
        self.isoma = []
        self.currents = {'ampa': [], 'nmda': []}
        self.all_releases = []
        self.all_release_events = []
        for nrep in xrange(iterations): # could do multiple runs.... 
            self.reset()
            self['v_pre'] = pre_cell.soma(0.5)._ref_v
            self['t'] = h._ref_t
            self['v_soma'] = pre_cell.soma(0.5)._ref_v
            self['relsite_xmtr'] = synapse.terminal.relsite._ref_XMTR[0]

            if isinstance(synapse.psd, GluPSD):
                # make a synapse monitor for each release zone
                self.all_nmda = []
                self.all_ampa = []
                for syn in synapses:
                    # collect all PSDs across all synapses
                    self.all_ampa.extend(syn.psd.ampa_psd)
                    self.all_nmda.extend(syn.psd.nmda_psd)
                    
                    # Record current through all PSDs individually
                    syn.psd.record('i', 'g', 'Open')
            
                #for k,p in enumerate(self.all_nmda):
                    #self['iNMDA%03d' % k] = p._ref_i
                    #self['opNMDA%03d' % k] = p._ref_Open
                #for k,p in enumerate(self.all_ampa):
                    #self['iAMPA%03d' % k] = p._ref_i
                    #self['opAMPA%03d' % k] = p._ref_Open
        
            elif isinstance(synapse.psd, GlyPSD):
                #  Record current through all PSDs individually
                for k,p in enumerate(self.all_psd):
                    self['iGLY%03d' % k] = p._ref_i
                    self['opGLY%03d' % k] = p._ref_Open
                
                psd = self.all_psd
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

            for i, s in enumerate(synapses):
                s.terminal.relsite.rseed = util.random.current_seed() + nrep
            self.custom_init()
            h.run()

            # add up psd current across all runs
            if isinstance(synapse.psd, GluPSD):
                iampa = np.zeros_like(synapse.psd.get_vector('ampa', 'i'))
                inmda = iampa.copy()
                for syn in self.synapses:
                    for i in range(syn.psd.n_psd):
                        iampa += syn.psd.get_vector('ampa', 'i', i)
                        inmda += syn.psd.get_vector('nmda', 'i', i)
                isoma = iampa + inmda
                self.currents['ampa'].append(iampa)
                self.currents['nmda'].append(inmda)
            elif isinstance(synapse.psd, GlyPSD):
                isoma = np.zeros_like(self['iGLY000'])
                for k in range(len(self.all_psd)):
                    isoma += self['iGLY%03d'%k]
            self.isoma.append(isoma)
            self.all_releases.append(self.release_timings())
            self.all_release_events.append(self.release_events())
            # if nrep > 0:
            #     assert (np.all(self.isoma[0] != self.isoma[1]))
#        print np.all(self.isoma[0] == self.isoma[1])
        
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

    def release_timings(self):
        """
        Return a list of arrays (one array per synapse) describing the timing 
        and latency of release events. 
        """
        data = []
        for j in range(0, len(self.synapses)):
            relsite = self.synapses[j].terminal.relsite
            nev = relsite.ev_index
            ev = np.empty(nev, dtype=[('time', float), ('latency', float)])
            ev['latency'] = np.array(relsite.EventLatencies)[:nev]
            ev['time'] = np.array(relsite.EventTime)[:nev]
            data.append(ev)
        return data
    
    def open_probability(self):
        """ 
        Analyze results and return a dict of values related to psd open 
        probability:
            
            nmda: (imax, opmax)
            ampa: (imax, opmax)
            gly:  (imax, opmax)
        """
        synapse = self.synapses[0]
        if isinstance(synapse.psd, GluPSD) and len(synapse.psd.nmda_psd) > 0:
            # find a psd with ampa and nmda currents
            nmImax = 0
            amImax = 0
            nmOmax = 0
            amOmax = 0
            #self.win.nextRow()
            for syn in self.synapses:
                for i in range(syn.psd.n_psd):
                    nm = np.abs(syn.psd.get_vector('nmda', 'i', i)).max()
                    am = np.abs(syn.psd.get_vector('ampa', 'i', i)).max()
                    opnm = np.abs(syn.psd.get_vector('nmda', 'Open', i)).max()
                    opam = np.abs(syn.psd.get_vector('ampa', 'Open', i)).max()
                    if nm > 1e-6 or am > 1e-6:
                        nmImax = nm
                        amImax = am
                        nmOmax = opnm
                        amOmax = opam
                        break
                if nmImax != 0:
                    break
            
            return {'nmda': OrderedDict([('Imax', nmImax), ('Omax', nmOmax)]), 
                    'ampa': OrderedDict([('Imax', amImax), ('Omax', amOmax)])}
        
        elif isinstance(synapse.psd, GlyPSD) and len(synapse.psd.all_psd) > 0:
            # find a psd with ampa and nmda currents
            glyImax = 0
            glyOmax = 0
            for i in range(len(self.all_psd)):
                imax = np.abs(self['iGLY%03d'%i]).max()
                omax = np.abs(self['opGLY%03d'%i]).max()
            
            return {'gly': (glyImax, glyOmax)}

    def analyze_events(self):
        events = []
        for run in range(len(self.isoma)):
            events.append(self.analyze_events_in_run(runno=run))
        return events
        
    def analyze_events_in_run(self, runno=0):
        """
        Analyze postsynaptic events for peak, latency, and shape.
        
        Todo: 
        - This currently analyzes cumulative currents; might be better to 
          analyze individual PSD currents
        - Measure decay time constant, rate of facilitation/depression,
          recovery.
        
        """
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
        events[:] = np.nan
        
        minLat = 0.0 # minimum latency for an event, in ms
        minStart = int(minLat / self.dt)  # first index relative to pulse to search for psc peak
        
        for i in range(stim['NP']):
            tstart = stim['delay'] + i * ipi # pulse start time 
            events['pulse time'][i] = tstart
            istart = int(tstart / self.dt)   # pulse start index
            tp[i] = tstart - stim['delay']
            iend = istart + pscpts
            #        print 'istart: %d iend: %d, len(isoma): %d\n' % (istart, iend, len(isoma))
            ipsc[i, :] = np.abs(self.isoma[runno][istart:iend])
            psc_pk = minStart + np.argmax(ipsc[i, minStart:-(extend_pts+1)]) # position of the peak
            
            #print 'i, pscpk, ipsc[i,pscpk]: ', i, psc_pk, ipsc[i, psc_pk]
            #       print 'minLat: %f   ipi+t_extend: %f, hdt: %f' % ((minLat, ipi+t_extend, self.dt))
            if psc_pk == minStart:
                continue
            pkval = ipsc[i, psc_pk]
            events['peak'][i] = pkval
            events['peak index'][i] = psc_pk
            
            # Find 20% and 80% crossing points to the left of the PSC peak
            pscmin = ipsc[i, :psc_pk].min()
            
            lat20 = find_crossing(ipsc[i], start=psc_pk, direction=-1,
                                  threshold=(pscmin + (pkval-pscmin) * 0.2)) * self.dt
            lat80 = find_crossing(ipsc[i], start=psc_pk, direction=-1,
                                  threshold=(pscmin + (pkval-pscmin) * 0.8)) * self.dt
            events['20% latency'][i] = lat20
            events['80% latency'][i] = lat80
            
            # Find 50% crossing points on either side of the PSC peak
            psc_50l = find_crossing(ipsc[i], start=psc_pk, direction=-1,
                                    threshold=(pscmin + (pkval-pscmin) * 0.5)) * self.dt
            psc_50r = find_crossing(ipsc[i], start=psc_pk, direction=1,
                                    threshold=(pscmin + (pkval-pscmin) * 0.5)) * self.dt
            
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

        return events

    def hide(self):
        if hasattr(self, 'win'):
            self.win.hide()

    def show(self, releasePlot=True, probabilityPlot=True, glyPlot=False, plotFocus='EPSC'):
        self.win = pg.GraphicsWindow()
        self.win.resize(1000, 1000)
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
        if 'gly' in oprob:
            glyImax, glyOPmax = oprob['gly']
            print 'Max GLYR Open Prob: %f' % (glyOPmax,)
            print 'Max GLYR I: %f' % (glyImax,)
        else:
            nmImax, nmOPmax = oprob['nmda'].values()
            amImax, amOPmax = oprob['ampa'].values()
            print 'Max NMDAR Open Prob: %f   AMPA Open Prob: %f' % (nmOPmax, amOPmax)
            print 'Max NMDAR I: %f   AMPA I: %f' % (nmImax, amImax)
            if nmImax + amImax != 0.0:
                print '   N/(N+A): %f\n' % (nmImax / (nmImax + amImax))
            else:
                print "   (no NMDA/AMPA current; release might have failed)"


        #
        # Plot pre/postsynaptic currents
        #
        t = self['t']

        p1 = self.win.addPlot(title=self.pre_cell.status['name'])
        p1.setLabels(left='V pre (mV)', bottom='Time (ms)')
        p1.plot(t, self['v_pre'])
        
        if plotFocus == 'EPSC':
            self.win.nextRow()
            p2 = self.win.addPlot(title=self.post_cell.status['name'])
            for i, isoma in enumerate(self.isoma):
                p2.plot(t, isoma, pen=(i, len(self.isoma)))
            p2.plot(t, np.mean(self.isoma, axis=0), pen=pg.mkPen('w', width=2))
            p2.setLabels(left='Total PSD current (nA)', bottom='Time (ms)')
            p2.setXLink(p1)
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
        
        eventno = 0
        self.win.nextRow()
        p3 = self.win.addPlot(labels={'left': '20%-80% Latency (ms)', 'bottom': 'Pulse Time (ms)'})
        p3.plot(events[eventno]['pulse time'], events[eventno]['20% latency'], pen=None, symbol='o')
        p3.plot(events[eventno]['pulse time'], events[eventno]['80% latency'], pen=None, symbol='t')
        p3.setXLink(p1)
        
        self.win.nextRow()
        p4 = self.win.addPlot(labels={'left': 'Half Width (ms)', 'bottom': 'Pulse Time (ms)'})
        p4.plot(events[eventno]['pulse time'], events[eventno]['half width'], pen=None, symbol='o')
        p4.setXLink(p1)
        
        self.win.nextRow()
        p5 = self.win.addPlot(labels={'left': 'Rise Time (ms)', 'bottom': 'Pulse Time (ms)'})
        p5.plot(events[eventno]['pulse time'], events[eventno]['rise time'], pen=None, symbol='o')
        p5.setXLink(p1)
        
        
        #
        # Print average values from events
        #
        nst = range(self.stim['NP'])
        analysisWindow = [nst[0:2], nst[-10:-1]]
        print analysisWindow
        print events[eventno]['rise time']
        RT_mean2080_early = np.nanmean(events[eventno]['rise time'][analysisWindow[0]])
        RT_mean2080_late = np.nanmean(events[eventno]['rise time'][analysisWindow[1]])
        Lat_mean20_early = np.nanmean(events[eventno]['20% latency'][analysisWindow[0]])
        Lat_mean20_late = np.nanmean(events[eventno]['20% latency'][analysisWindow[1]])
        HW_mean_early = np.nanmean(events[eventno]['half width'][analysisWindow[0]])
        HW_mean_late = np.nanmean(events[eventno]['half width'][analysisWindow[1]])
        print "\n--------------"
        print "Means:"
        print "--------------"
        #print RT_mean2080_early
        #print Lat_mean20_early
        #print HW_mean_early
        print 'Early:   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_mean2080_early, Lat_mean20_early,
                                                                                HW_mean_early)
        print 'Late :   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_mean2080_late, Lat_mean20_late,
                                                                                HW_mean_late)
        RT_std2080_early = np.nanstd(events[eventno]['rise time'][analysisWindow[0]])
        RT_std2080_late = np.nanstd(events[eventno]['rise time'][analysisWindow[1]])
        Lat_std20_early = np.nanstd(events[eventno]['20% latency'][analysisWindow[0]])
        Lat_std20_late = np.nanstd(events[eventno]['20% latency'][analysisWindow[1]])
        HW_std_early = np.nanstd(events[eventno]['half width'][analysisWindow[0]])
        HW_std_late = np.nanstd(events[eventno]['half width'][analysisWindow[1]])
        print "\n--------------"
        print "Standard Deviations:"
        print "--------------"
        print 'Early:   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_std2080_early, Lat_std20_early,
                                                                                HW_std_early)
        print 'Late :   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_std2080_late, Lat_std20_late,
                                                                                HW_std_late)
        print "-----------------"


        #
        # Plot release event distributions over time
        #
        if releasePlot:
            self.win.nextRow()
            p6 = self.win.addPlot(labels={'left': 'Release latency', 'bottom': 'Time (ms)'})
            p6.setXLink(p1)
            p7 = self.win.addPlot(labels={'left': 'Release latency', 'bottom': 'Num. Releases'})
            p7.setYLink(p6)
            self.win.ci.layout.setColumnFixedWidth(1, 200)
            
            events = self.all_releases
            all_latencies = []
            for i, trial in enumerate(events):
                for syn in trial:
                    p6.plot(syn['time'], syn['latency'], pen=None, symbolBrush=(i, len(events)), 
                        symbolPen=(i, len(events)), symbolSize=4, symbol='o')
                    all_latencies.append(syn['latency'])
            all_latencies = np.concatenate(all_latencies)
            (hist, binedges) = np.histogram(all_latencies)
            curve = p7.plot(binedges, hist, stepMode=True, fillBrush=(100, 100, 255, 150), fillLevel=0)
            curve.rotate(-90)
            curve.scale(-1, 1)

        if probabilityPlot:
            self.win.nextRow()
            p8 = self.win.addPlot(labels={'left': 'Release Prob', 'bottom': 'Time (ms)'})
            p8.setXLink(p1)
            events = self.all_releases

            # for i, trial in enumerate(events):
            #     for syn in trial:
            #         p6.plot(syn['time'], syn['latency'], pen=None, symbolBrush=(i, len(events)),
            #             symbolPen=(i, len(events)), symbolSize=4, symbol='o')
            #         all_latencies.append(syn['latency'])
            
        #
        # Plot GlyR state variables
        #
        # if glyPlot:
       #      i = 0
       #      if synapse.psd.psdType == 'glyslow':
       #          mpl.figure(2)
       #          for var in ['C0', 'C1', 'C2', 'O1', 'O1', 'D1', 'Open']:
       #              mpl.subplot(nstate, 1, i + 1)
       #              mpl.plot(t, self[var])
       #              mpl.ylabel(var)
       #              i = i + 1
       #      if synapse.psd.psdType == 'glyfast':
       #          mpl.figure(2)
       #          for var in ['C0', 'C1', 'C2', 'C3', 'O1', 'O2', 'Open']:
       #              mpl.subplot(7, 1, i + 1)
       #              mpl.plot(t, self[var])
       #              mpl.ylabel(var)
       #              i = i + 1
       #      mpl.draw()
        #mpl.show()
