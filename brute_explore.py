#!/usr/bin/python

"""
Brute force parameter space exploration

Two steps:
1. Rin, taum and Vm exploration
2. FI curve exploration (once resting parameters are set)
"""
import os
import sys
import numpy as np

from neuron import h

import cnmodel
import cnmodel.cells as cells
from cnmodel.protocols import IVCurve, VCCurve
from cnmodel.util import nstomho
from cnmodel.util.stim import make_pulse
import pyqtgraph as pg


path = os.path.dirname(cnmodel.__file__)
h.load_file("stdrun.hoc")
h.load_file(os.path.join(path, "custom_init.hoc"))

def spike_times(Vm, dt, threshold=0.):
    """
    Return an array of spike times for each trace.
    Works only on ONE trace in this version...
    
    Parameters
    ==========
    threshold: float (default: None)
        Optional threshold at which to detect spikes. By 
        default, this queries cell.spike_threshold.
    
    Returns
    =======
    list of spike times.
    
    """
    Vm = np.array(Vm)
#    spikes = []
    #dvdt = np.diff(Vm[i]) / self.dt
    #mask = (dvdt > 40).astype(int)
    mask = (Vm > threshold).astype(int)
    indexes = np.argwhere(np.diff(mask) == 1) + 1
    times = indexes.astype(float) * dt
#    spikes.append(times)
    return times

def find_minimumFI(target, slmap):
    """
    Find the minimum distance between points in the parmap and the target
    algorithm is simple:
    targetfr - fr (minimum)
    
    Parameters
    ----------
    target: dict
        {'fr': vmval}
    
    """
    # weights
    wFR = 1.0
    diffs = np.zeros(len(slmap))
    for i, s in enumerate(slmap):
        fr = s['spikes']
        frd = np.fabs(fr-target['fr'])
        diffs[i] = np.sqrt((wFR*frd)**2) #  + wRin*rind**2 + wtau*taud**2)
    best = np.argmin(diffs)
    return slmap[best]['pars'], slmap[best]['spikes']
    
def find_minimum(target, slmap):
    """
    Find the minimum distance between points in the parmap and the target
    algorithm is simple:
    Vtarget - vm (minimum)
    1.0 - Rintarget/Rin (value closest to 0)
    1.0 - tautarget/tau (value closest to 0)
    
    Parameters
    ----------
    target: dict
        {'Vm': vmval, 'Rin': rinval, 'tau': tval}
    
    """
    # weights
    wVm = 1.0
    wRin = 1.0
    wtau = 1.0
    diffs = np.zeros(len(slmap))
    for i, s in enumerate(slmap):
        vm = s['vals']['v']
        rin = s['vals']['Rin']
        tau = s['vals']['tau']
        vmd = np.fabs(vm-target['vm'])
        rind = 1.0 - (rin/target['Rin'])
        taud = 1.0 - (tau/target['tau'])
        diffs[i] = np.sqrt(wVm*vmd**2 + wRin*rind**2 + wtau*taud**2)
    print diffs
    best = np.argmin(diffs)
    return slmap[best]['pars'], slmap[best]['vals']
 
def brute_search_rmtauvm(cell):
    print 'setup'
    erevs = np.arange(-90., -50., 1.)
    leaks = np.arange(0., 10.1, 0.1)
    ihgbar = np.arange(0.5, 6.1, 0.05)
    nsomas = len(erevs)
    nleaks = len(leaks)
    nih = len(ihgbar)
    slmap = []
    cell.vrange=[-100., 100.]
    for erev in erevs:
        cell.soma().leak.erev = erev
        #cell.set_soma_size_from_Cm(somasize)
        for leakbar in leaks:
            cell.soma().leak.gbar = nstomho(leakbar, cell.somaarea)
            for ihbar in ihgbar:
                cell.soma().ihvcn.gbar = nstomho(ihbar, cell.somaarea)
                cell.vm0 = None
                cell.get_mechs(cell.soma)
                cell.vrange=[-100., 100.]
                h.celsius = 35
                h.batch_run(20., 0.025, "cellinit.dat")
                
#                cell.cell_initialize(vrange=cell.vrange)
                resting = cell.compute_rmrintau(auto_initialize=True, vrange=cell.vrange)
                sys.stdout.write("\r" + 'erev: %5.1f leakbar: %6.3f, ihbar: %6.3f  %3d' %
                     (erev, cell.soma().leak.gbar*1000., cell.soma().ihvcn.gbar*1000., len(slmap)))
                sys.stdout.write(' Vm: %5.1f Rin: %6.1f, tau: %6.1f  %3d' %
                     (resting['v'], resting['Rin'], resting['tau'], len(slmap)))
                     
                if  ((8. < resting['tau'] < 12.) and
                    (-73. < resting['v'] < -67.) and
                    (120. < resting['Rin'] < 170.)):
                    slmap.append({'vals': resting, 'pars': {'soma': somasize, 'ihbar': ihbar, 'leakbar': leakbar}})
                    
                    # print ('somasize: %.1f leakbar: %.3f, ihbar: %.3f v: %.2f, Rin: %.1f  tau: %.1f' %
                    # (somasize, cell.soma().leak.gbar*1000., cell.soma().ihvcn.gbar*1000., resting['v'], resting['Rin'], resting['tau']))
                #print '   L, h: ', cell.soma().cm, cell.soma().area(), cell.soma().diam, cell.soma().sec.L
            # except:
            #     pass
    print '\n----------------'
    return slmap

def tune_fi(cell, iinj, Vplot=None):
    """
    Given basic parameters, adjust nachans ka and kht for firing rate at specific current in 100 msec pulse.
    """
    nabars = np.arange(800., 1600., 20.)
    khtbars = np.arange(100., 400., 20.)
    kabars = np.arange(65., 250., 250.)
    dt = 0.025
    tend = 120.
    stim = {
        'NP': 1,
        'delay': 5.,
        'dur': 100.,
        'amp': 1.0,
        'dt': dt,
        }
    istim = h.iStim(0.5, sec=cell.soma)
    istim.delay = 5.
    istim.dur = 1e9 # these actually do not matter...
    istim.iMax = 0.0
    stim['amp'] = iinj
    nostim = {
        'NP': 1,
        'delay': 5.,
        'dur': 100.,
        'amp': 0.0,
        'dt': dt,
        }    
    
    (secmd, maxt, tstims) = make_pulse(stim)
    (nostim, maxt, tstims) = make_pulse(nostim)
    vecs = {'v_soma': h.Vector(), 'i_inj': h.Vector(), 'time': h.Vector()}
    vecs['i_stim'] = h.Vector(secmd)
    vecs['nostim'] = h.Vector(nostim)
    
    slmap = []
    rn = 0
    nlst = 0
    for ina, nabar in enumerate(nabars):
        cell.adjust_na_chans(cell.soma, gbar=nabar)
        cxs = len(nabars)
        colors = [(ic, cxs*3./2.) for ic in range(cxs)]
        for j, khtbar in enumerate(khtbars):
            cell.soma().kht.gbar = nstomho(khtbar, cell.somaarea)
            for kabar in kabars:
                cell.soma().ka.gbar = nstomho(kabar, cell.somaarea)
                vecs['v_soma'].record(cell.soma(0.5)._ref_v)
                vecs['i_inj'].record(istim._ref_i)
                vecs['time'].record(h._ref_t)
                
                vecs['nostim'].play(istim._ref_i, h.dt, 0, sec=cell.soma)
                h.batch_run(100., 0.025, "cellinit.dat")                

                cell.cell_initialize()  # initialize the cell to it's rmp
                
                # connect current command vector
                vecs['i_stim'].play(istim._ref_i, h.dt, 0, sec=cell.soma)
                # GO
                h.dt = dt
                h.celsius = 35
                h.tstop = tend
                cell.vm0 = None
                r = cell.compute_rmrintau(auto_initialize=False)
                Rin, tau, v = r['Rin'], r['tau'], r['v']
#                print '    *** Rin: %9.0f  tau: %9.1f   v: %6.1f' % (Rin, tau, v)

                cell.cell_initialize(vrange=cell.vrange)  # initialize the cell to it's rmp
                while h.t < h.tstop:
                    h.fadvance()
                st = spike_times(np.array(vecs['v_soma']), dt, threshold=-20.)
#                sys.stdout.write("%d" % rn)
                rn += 1
                if Vplot is not None:
                    Vplot.plot(vecs['time'], vecs['v_soma'], pen=colors[ina])
                if len(st) > 2:
                    sys.stdout.write("\r" + 'spks: %d  gna: %.1f ght: %.3f, ga: %.3f  n=%d' %
                         (len(st), cell.soma().na.gbar*1000., cell.soma().kht.gbar*1000., cell.soma().ka.gbar*1000., len(slmap)))
                
                    slmap.append({'spikes': len(st), 'pars': {'na': nabar, 'kht': khtbar, 'ka': kabar}})

    print('space explored: %d found %d\n' % (rn, len(slmap)))
    return slmap


if __name__ == '__main__':
    target = {'vm': -72., 'Rin': 142., 'tau': 10.}
    cellobj = cells.Tuberculoventral
    cell = cellobj.create(debug=False, ttx=False, modelType='TVmouse', species='mouse', morphology='cnmodel/morphology/tv_stick.hoc', decorator=True)
#    cell = cellobj.create(debug=False, ttx=False, modelType='TVmouse', species='mouse')
    cell.set_d_lambda(freq=2000) 
    # init the cell (should do this inside loops...)
    h.topology()
    if sys.argv[1] == 'rm':
        slmap = brute_search_rmtauvm(cell)
        pars, vals = find_minimum(target, slmap)
        print pars
        print vals

    if sys.argv[1] == 'fi':
        app = pg.mkQApp()
        win = pg.GraphicsWindow('hi there')
        win.resize(1000, 800)
        Vplot = win.addPlot(labels={'left': 'Vm (mV)', 'bottom': 'Time (ms)'})
        slmap = tune_fi(cell, 0.600, Vplot)
        pg.show()
        if sys.flags.interactive == 0:
            pg.QtGui.QApplication.exec_()
        pars, vals = find_minimumFI({'fr': 40.}, slmap)
        print pars
        print vals


    