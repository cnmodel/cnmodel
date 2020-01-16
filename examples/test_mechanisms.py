"""
test_mechanisms.py

This program displays the results of inserting NEURON mechaanisms from .mod files 
into a point cell to voltage steps in voltage clamp.
This code is primarily for visual verification of model function.

Usage: python test_mechanisms.py <mechname>

Available mechanisms::

   CaPCalyx                bkpkj                 hcno               hcnobo                 hpkj
      ihpyr          ihsgcApical     ihsgcBasalMiddle                ihvcn                jsrna
         ka                 kcnq                kdpyr                  kht                  kif
        kis                  klt                 kpkj                kpkj2             kpkjslow
      kpksk                 leak                lkpkj                   na                naRsg
       nacn             nacncoop                  nap                napyr                nav11
      nacsh

Note: only modfiles that implement voltage-dependent ion channel models make sense to run
with this routine. the list "nottestablemechs" in the file defines mechanisms provided
with cnmodel that cannot be run with this program.

"""
import sys
from neuron import h
from neuron import nrn

import gc
import numpy as np
#import scipy as sp
import cnmodel.util
import pyqtgraph as pg
import pyqtgraph.exporters
from pyqtgraph.Qt import QtCore, QtGui
import cnmodel.util.pynrnutilities as Util

nottestablemechs = ['cadyn', 'ca_ion', 'cadiff', 'cadifpmp', 'Mechanism',
                'capmp', 'capump', 'cl_ion', 'extracellular', 'fastpas',
                'k_ion', 'KIR', 'hh', 'na_ion', 'narsg', 'pas', 'cap']  # cap uses "pcabar"

class ChannelKinetics():
    def __init__(self, args, export=False):
        modfile = []
        if isinstance(args, list):
            for arg in args:
                modfile.append(arg) # must be string, not list...
        else:
            modfile.append(args) # 'CaPCalyx'
        print('modfile: ', modfile)
        colors = ['w', 'r', 'g', 'y', 'c', 'm', 'b']
        if len(modfile) > len(colors):
            print('Too many modfiles... keep it simple!')
            exit()
        # if isinstance(args, list) and len(args) > 1:
        #     modfile2 = args[1]
        doKinetics = False
        self.app = pg.mkQApp()
        self.win = pg.GraphicsWindow()
        self.win.setWindowTitle('VC Plots')
        self.win.resize(900,600)
        # cw = QtGui.QWidget()
        # self.win.setCentralWidget(cw)
        # self.gridLayout = QtGui.QGridLayout()
        # cw.setLayout(self.gridLayout)
        # self.gridLayout.setContentsMargins(9, 9, 4, 4)
        # self.gridLayout.setSpacing(1)
        self.p1 = self.win.addPlot(title="I (VC)")
        # self.gridLayout.addWidget(self.p1, 0, 0, 1, 1)
        self.p2 = self.win.addPlot(title="I_ss, I_max")
        # self.gridLayout.addWidget(self.p2, 0, 1, 1, 1)
        self.p2mh = self.win.addPlot(title="act, inact")
        
        self.win.nextRow()
        self.p3 = self.win.addPlot(title="V command")
        # self.gridLayout.addWidget(self.p3, 1, 0, 1, 1)
        self.p5 = self.win.addPlot(title="I_min")
        # self.gridLayout.addWidget(self.p5, 1, 1, 1, 1)
        self.p6 = self.win.addPlot(title='Inact')
        self.win.show()

        QtGui.QApplication.processEvents()
        #
        # self.tdur is a table of durations for the pulse and post-pulse for each channel type (best to highlight features
        # on appropriate time scales)
        #
        self.tdur = {'CaPCalyx': [20., 10.], 
                     'nav11': [10., 5.], 'jsrna': [10., 5.], 'ichanWT2005': [10., 5.], 
                     'nacn': [10., 5.], 'nacncoop': [10., 5.], 'nabu': [10., 5.], 'nacsh': [10., 5.],
                     'kht':[200., 20.], 'klt': [200., 20.], 'ka': [25., 5.],
                     'hcno': [1000., 200.], 'ih': [1000., 200.], 'ihvcn': [1000., 200.],'hcnobo': [1000., 200.],
                     'ihsgcBasalMiddle': [1000., 200.], 'ihsgcApical': [1000., 200.], 
                     'kif': [100., 100.], 'kis': [100., 10.], 'napyr': [10, 5.], 'ihpyr': [1000., 200.],
                     'kdpyr': [200., 20.], 'kcnq': [200, 20], 'nap': [200., 100.],
                     }
        for i, mfile in enumerate(modfile):
            self.run(modfile=mfile, color=colors[i], export=export)
        s=''
        #self.win.setWindowTitle('VC Plots: ' + [s+sn+'; ' for sn in modfile])
        gc.collect()

        if doKinetics:
            self.win2 = pg.GraphicsWindow(title='KineticPlots')
            self.win2.resize(800, 600)
            self.kp1 = self.win.addPlot(title="htau")
            self.computeKinetics('nav11')
            
        #self.win.show()

    def run(self, modfile='CaPCalyx', color='r', export=False):
        if isinstance(modfile, list):
            modfile = modfile[0]

        if modfile in self.tdur:
            tstep = self.tdur[modfile]
        else:
            tstep = [200., 50.]
        
        tdelay = 5.0
        Channel = cnmodel.util.Mechanism(modfile)
        leak = cnmodel.util.Mechanism('leak')
        if modfile == 'nacsh':
            Channel.set_parameters({'gbar': 6e3})  # weird units from Kole paper
        else:
            Channel.set_parameters({'gbar': 1})
        leak.set_parameters({'gbar': 1e-12})

        self.soma = cnmodel.util.Section(L=10, diam=10, mechanisms=[Channel, leak])
        if modfile == 'bkpjk':
            ca_init = 100e-6
            self.soma().cai = ca_init
        else:
            ca_init = 70e-6
        if modfile == 'nacncoop':
            self.soma().nacncoop.p = 0.1
            self.soma().nacncoop.KJ = 1000.
        Vr = 0.
        if modfile.startswith('na') or modfile.startswith('jsr') or modfile.startswith('ichanWT'):
            Vr = 50.
        elif modfile.startswith('k'):
            Vr = -84.
        elif modfile.startswith('ih') or modfile.startswith('hcn'):
            Vr = -43.
        elif modfile.startswith('Ca') or modfile.startswith('ca'):
            Vr = +100.
        h.celsius = 37. # set the temperature.
        self.vec={}
        for var in ['time', 'V', 'IChan', 'Vcmd']:
            self.vec[var] = h.Vector()

        h.dt = 0.025
        v_init = -65.

        clampV = v_init
        self.vcPost = h.SEClamp(0.5, sec=self.soma)
        self.vcPost.dur1 = tdelay
        self.vcPost.amp1 = clampV
        self.vcPost.dur2 = tstep[0]
        self.vcPost.amp2 = clampV-0.0 # just a tiny step to keep the system honest
        self.vcPost.dur3 = tstep[1]
        self.vcPost.amp3 = clampV
        self.vcPost.rs = 1e-9
        print("soma: ", self.soma, end=' ') 
        # print(dir(self.vcPost))
        # print(' vcpost sec: ', self.vcPost.Section())

        if modfile[0:2] == 'ih':
            stimamp = np.linspace(-140, -40, num=21, endpoint=True)
        else:
            stimamp = np.linspace(-100, 40, num=35, endpoint=True)
        self.ivss = np.zeros((2, stimamp.shape[0]))
        self.ivmin = np.zeros((2, stimamp.shape[0]))
        self.ivmax = np.zeros((2, stimamp.shape[0]))
        print(('I range = %6.1f-%6.1f, T = %4.1f' % (np.min(stimamp), np.max(stimamp), h.celsius)))

        for i, V in enumerate(stimamp):
            self.vcPost.dur1 = tdelay
            self.vcPost.amp1 = clampV
            self.vcPost.dur2 = tstep[0]
            self.vcPost.dur3 = tstep[1]
            self.vcPost.amp3 = clampV
            self.vcPost.amp2 = V
            self.vec['IChan'].record(self.vcPost._ref_i, sec=self.soma)
            self.vec['V'].record(self.soma()._ref_v, sec=self.soma)
            self.vec['time'].record(h._ref_t)
#            print 
            h.tstop = self.vcPost.dur1+self.vcPost.dur2+self.vcPost.dur3
            h.finitialize(v_init)
            h.run()
            self.t = np.array(self.vec['time'])
            self.ichan = np.array(self.vec['IChan'])
            self.v = np.array(self.vec['V'])
            self.p1.plot(self.t, self.ichan, pen=pg.mkPen((i, len(stimamp)*1.5)))
            self.p3.plot(self.t, self.v, pen=pg.mkPen((i, len(stimamp)*1.5)))
            (self.ivss[1, i], r2) = Util.measure('mean', self.t, self.ichan, tdelay+tstep[0]-10., tdelay+tstep[0])
            (self.ivmin[1, i], r2) = Util.measure('min', self.t, self.ichan, tdelay+0.1, tdelay+tstep[0]/5.0)
            (self.ivmax[1, i], r2) = Util.measure('max', self.t, self.ichan, tdelay+0.1, tdelay+tstep[0]/5.0)
            self.ivss[0,i] = V
            self.ivmin[0,i] = V
            self.ivmax[0,i] = V
        self.p2.plot(self.ivss[0,:], np.array(self.ivss[1,:]), pen=pg.mkPen(color), symbol='o', symbolSize=4.0, )
        self.p2.plot(self.ivmax[0,:], np.array(self.ivmax[1,:]), pen=pg.mkPen(color), symbol='t', symbolSize=4.0, )
        self.p5.plot(self.ivmin[0,:], np.array(self.ivmin[1,:]), pen=pg.mkPen(color), symbol='s', symbolSize=4.0, )
        # activation
        color_act = 'y'
        color_inact = 'c'
        tdelay = 200.
        tstep[0] = 100.
        tstep[1] = 50.
        for i, V in enumerate(stimamp):
            self.vcPost.dur1 = tdelay
            self.vcPost.amp1 = -120.
            self.vcPost.dur2 = tstep[0]
            self.vcPost.amp2 = V
            self.vcPost.dur3 = tstep[1]
            self.vcPost.amp3 = clampV
            self.vec['IChan'].record(self.vcPost._ref_i, sec=self.soma)
            self.vec['V'].record(self.soma()._ref_v, sec=self.soma)
            self.vec['time'].record(h._ref_t)
#            print
            h.tstop = self.vcPost.dur1+self.vcPost.dur2+self.vcPost.dur3
            h.finitialize(v_init)
            h.run()
            self.t = np.array(self.vec['time'])
            self.ichan = np.array(self.vec['IChan'])
            self.v = np.array(self.vec['V'])
            # self.p1.plot(self.t, self.ichan, pen=pg.mkPen((i, len(stimamp)*1.5)))
            # self.p3.plot(self.t, self.v, pen=pg.mkPen((i, len(stimamp)*1.5)))
            (self.ivss[1, i], r2) = Util.measure('mean', self.t, self.ichan, tdelay+tstep[0]-10., tdelay+tstep[0])
            (self.ivmin[1, i], r2) = Util.measure('min', self.t, self.ichan, tdelay+0.1, tdelay+tstep[0]/5.0)
            (self.ivmax[1, i], r2) = Util.measure('max', self.t, self.ichan, tdelay+0.1, tdelay+tstep[0]/5.0)
            self.ivss[0,i] = V
            self.ivmin[0,i] = V
            self.ivmax[0,i] = V
        if Vr > 0:
            fm = self.ivmin
        else:
            fm = self.ivmax
        # self.p2mh.plot(fm[0,:], self.ivmax[1,:], symbol='t', symbolSize=4.0, pen=pg.mkPen(color_act))
        G = fm[1,:]/(fm[0,:]-Vr)
        Po = G/np.max(G)
        self.p2mh.plot(fm[0,:], Po, symbol='s', symbolSize=4.0, pen=pg.mkPen(color_act))

        # inactivation
        tdelay = 200.
        tstep[0] = 100.
        tstep[1] = 50.
        vstep = 0.
        for i, V in enumerate(stimamp):
            v_init = V
            self.vcPost.dur1 = tdelay
            self.vcPost.amp1 = V
            self.vcPost.dur2 = tstep[0]
            self.vcPost.amp2 = 0.
            self.vcPost.dur3 = tstep[1]
            self.vcPost.amp3 = clampV
            self.vec['IChan'].record(self.vcPost._ref_i, sec=self.soma)
            self.vec['V'].record(self.soma()._ref_v, sec=self.soma)
            self.vec['time'].record(h._ref_t)
 #            print
            h.tstop = self.vcPost.dur1+self.vcPost.dur2+self.vcPost.dur3
            h.finitialize(v_init)
            h.run()
            self.t = np.array(self.vec['time'])
            self.ichan = np.array(self.vec['IChan'])
            self.v = np.array(self.vec['V'])
            self.p6.plot(self.t[2:], self.ichan[2:], pen=pg.mkPen((i, len(stimamp)*1.5)))
            # self.p6.plot(self.t, self.v, pen=pg.mkPen((i, len(stimamp)*1.5)))
            (self.ivss[1, i], r2) = Util.measure('mean', self.t, self.ichan, tdelay+tstep[0]-10., tdelay+tstep[0])
            (self.ivmin[1, i], r2) = Util.measure('min', self.t, self.ichan, tdelay+0.1, tdelay+tstep[0]/5.0)
            (self.ivmax[1, i], r2) = Util.measure('max', self.t, self.ichan, tdelay+0.1, tdelay+tstep[0]/5.0)
            self.ivss[0,i] = V
            self.ivmin[0,i] = V
            self.ivmax[0,i] = V
        # self.p2mh.plot(self.ivmax[0,:], np.array(self.ivmax[1,:]), symbol='t', symbolSize=4.0, pen=pg.mkPen(color_inact))
        if Vr > 0:
            fm = self.ivmin
        else:
            fm = self.ivmax

        G = fm[1,:]/(vstep-Vr)
        Po = np.array(G/np.max(G))
        self.p2mh.plot(fm[0,:], Po, symbol='s', symbolSize=4.0, pen=pg.mkPen(color_inact))

                
        print(export)
        if export:
            exporter = pg.exporters.MatplotlibExporter(self.p1)
            print('exporting: ' + '%s_traces.svg' % modfile)
            exporter.export(fileName='%s_traces.pdf' % modfile)
            exporter = pg.exporters.MatplotlibExporter(self.p3)
            exporter.export('%s_command.pdf' % modfile)
            exporter = pg.exporters.MatplotlibExporter(self.p2)
            exporter.export('%s_IV.pdf' % modfile)
            
    def computeKinetics(self, ch):
        pass

def getmechs():
    mechs = []
    for n in dir(nrn):
        o = getattr(nrn, n)
        if str(o) == str(nrn.Mechanism):
            mechs.append(n)
    return mechs

def main():
    mechs = getmechs()
    if len(sys.argv) < 2:

        print("\n\nUsage: python test_mechanisms.py <mechname>")
        print("  Available mechanisms:")

        linelen = 0
        for i, n in enumerate(mechs):
            if n in nottestablemechs: # 'Mechanism':
                continue
            print("%20s" % n, end=' ')
            linelen += 20
            if linelen > 80:
                print("")

                linelen = 0
        sys.exit(1)
    
    if sys.argv[1:] in nottestablemechs:
        exit()
    export = False
    if len(sys.argv) > 2:
        if sys.argv[2] == 'export':
            export = True
    
    if sys.argv[1] == 'all':
        for n in mechs:
            if n in nottestablemechs:
                print(('Skipping %s' % n))
                continue
            else:
                ck = ChannelKinetics(n)
    else:
        ck = ChannelKinetics(sys.argv[1], export=export)
     
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()

if __name__ == "__main__":
    main()
