"""
Test synaptic connections between two different cell types. 

Usage:  python compare_simple_multisynapses.py <pre_celltype> <post_celltype>

This script:

1. Creates single pre- and postsynaptic cells
2. Creates a single synaptic terminal between the two cells, using the multisite synapse method.
3. Stimulates the presynaptic cell by current injection.
4. Records and analyzes the resulting post-synaptic events.
5. Repeats 3, 4 100 times to get an average postsynaptic event.
6. stores the resulting waveform in a pandas database

This is used mainly to check that the strength, kinetics, and dynamics of 
each synapse type is working as expected. 
"""

import numpy as np
import pyqtgraph as pg
import matplotlib.pyplot as mpl
import cnmodel.util.PlotHelpers as PH
from cnmodel.protocols import SynapseTest
from cnmodel import cells
from cnmodel.synapses import Synapse
import pickle
import lmfit

import sys

args = sys.argv[1:]

convergence = {
    'sgc': {'bushy': 1, 'tstellate': 1, 'dstellate': 1,  'octopus': 1, 'pyramidal': 1, 'tuberculoventral': 1},
    'dstellate': {'bushy': 1, 'tstellate': 1, 'dstellate': 1, 'pyramidal': 1, 'tuberculoventral': 1},
    'tuberculoventral': {'bushy': 1, 'tstellate': 1, 'pyramidal': 1, 'tuberculoventral': 1},
    }


class Exp2SynFitting:
    """
    Fit waveform against Exp2SYN function (used by Neuron)
    THe function is (from the documentation):
    i = G * (v - e)      i(nanoamps), g(micromhos);
    G = weight * factor * (exp(-t/tau2) - exp(-t/tau1))
    WHere factor is evaluated when initializing(in neuron; see exp2syn.mod source) as:
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)  (log10)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
    
    Parameters
    ----------
    initpars : dict
        dict of initial parameters. For example: {'tau1': 0.1, 
        'tau2': 0.3, 'weight': 0.1, 'delay' : 0.0, 'erev': -80.} (erev in mV)
    bounds : dict
        dictionary of bounds for each parameter, with a list of lower and upper values.
    
    """
    def __init__(self, initpars=None, bounds=None, functype='exp2syn'):
        self.fitpars = lmfit.Parameters()
        if functype == 'exp2syn':  # future handle other functions like alpha
            # (Name,  Value,  Vary,   Min,  Max,  Expr)
            self.fitpars.add_many(('tau1', initpars['tau1'], True, 0.05, 25., None),
                                  ('tau2', initpars['tau2'], True, 0.1, 50., None),
                                  ('weight', initpars['weight'], True, 1e-6, 1, None),
                                  ('erev', initpars['erev'], False), # do not adjust!
                                  ('v', initpars['v'], False),
                                  ('delay', initpars['delay'], True, 0., 5., None))
            self.func = self.exp2syn_err
        else:
            raise ValueError

    def fit(self, x, y, p, verbose=False):
        
        kws={'maxfev': 5000}
        # print('p: ', p)
        mim = lmfit.minimize(self.func, p, method='least_squares', args=(x, y)) #, kws=kws)
        if verbose:
            lmfit.printfuncs.report_fit(mim.params)
        fitpars = mim.params
        return fitpars

    #@staticmethod
    def exp2syn(self, x, tau1, tau2, weight, erev, v, delay):
        if tau1/tau2 > 1.0:
            tau1 = 0.999*tau2
        tp = (tau1*tau2)/(tau2 - tau1) * np.log10(tau2/tau1)
        factor = -np.exp(-tp/tau1) + np.exp(-tp/tau2)
        factor = 1.0/factor
        G = weight * factor * (np.exp(-(x-delay)/tau2) - np.exp(-(x-delay)/tau1))
        G[x-delay < 0] = 0.
        i = G * (v - erev)  #      i(nanoamps), g(micromhos);
        return i
    
    def exp2syn_err(self, p, x, y):
        return np.fabs(y-self.exp2syn(x, **dict([(k,p.value) for k,p in p.items()])))

def testexp():
    F = Exp2SynFitting(initpars={'tau1': 0.2, 'tau2': 0.4, 'weight': 0.1, 'erev': -80., 'v': -60., 'delay': 0})
    t = np.arange(0, 10., 0.01)
    p = F.fitpars
    target = F.exp2syn(t, 0.1, 0.3, 0.01, -80, -60, 0)
    # print(F.fitpars)
    pars = F.fit(t, target, F.fitpars)
    print('fitting results: ', pars)


def runtest(args, synapsetype='multisite'):
    if len(args) < 2:
        print("Usage:  python compare_simple_multisynapses.py <pre_celltype> <post_celltype>")
        print("   Supported cell types: sgc, bushy, tstellate, dstellate, octopus, tuberculoventral, pyramidal, cartwheel")
        sys.exit(1)

    celltypes = args[:2]    
    
    c = []
    for cellType in celltypes:
        if cellType == 'sgc':
            cell = cells.SGC.create()
        elif cellType == 'tstellate':
            cell = cells.TStellate.create(debug=True, ttx=False)
        elif cellType == 'dstellate': # Type I-II Rothman model, similiar excitability (Xie/Manis, unpublished)
            cell = cells.DStellate.create(model='RM03', debug=True, ttx=False)
        # elif cellType == 'dstellate_eager': # From Eager et al.
        #     cell = cells.DStellate.create(model='Eager', debug=True, ttx=False)
        elif cellType == 'bushy':
            cell = cells.Bushy.create(debug=True, ttx=True)
        elif cellType == 'octopus':
            cell = cells.Octopus.create(debug=True, ttx=True)
        elif cellType == 'tuberculoventral':
            cell = cells.Tuberculoventral.create(debug=True, ttx=True)
        elif cellType == 'pyramidal':
            cell = cells.Pyramidal.create(debug=True, ttx=True)
        elif cellType == 'cartwheel':
            cell = cells.Cartwheel.create(debug=True, ttx=True)
        else:
            raise ValueError("Unknown cell type '%s'" % cellType)
        c.append(cell)

    preCell, postCell = c
    
    print('celltypes: ', celltypes)
    nTerminals = convergence.get(args[0], {}).get(args[1], None)
    if nTerminals is None:
        nTerminals = 1
        print("Warning: Unknown convergence for %s => %s, assuming %d" % (args[0], args[1], nTerminals))

    if args[1:3] == ['sgc', 'bushy']:
        niter = 5
    else:
        niter = 20
    assert(synapsetype in ['simple', 'multisite'])
    st = SynapseTest()
    dt = 0.025
    stim = {
        'NP': 1,
        'Sfreq': 100.0,
        'delay': 0.0,
        'dur': 0.5,
        'amp': 10.0,
        'PT': 0.0,
        'dt': dt,
    }
    st.run(preCell.soma, postCell.soma, nTerminals, dt=dt, vclamp=-65., iterations=niter, synapsetype=synapsetype,
        tstop = 50., stim_params=stim)
    # st.show_result()
    # st.plots['VPre'].setYRange(-70., 10.)
    # st.plots['EPSC'].setYRange(-2.0, 0.5)
    # st.plots['latency2080'].setYRange(0., 1.0)
    # st.plots['halfwidth'].setYRange(0., 1.0)
    # st.plots['RT'].setYRange(0., 0.2)
    # st.plots['latency'].setYRange(0., 1.0)
    # st.plots['latency_distribution'].setYRange(0., 1.0)
    return st  # need to keep st alive in memory

def fit_one(st, stk):
    if stk[0] in ['sgc', 'tstellate', 'granule']:
        erev = 0.
    else:
        erev = -80.
    F = Exp2SynFitting(initpars={'tau1': 0.2, 'tau2': 0.4, 'weight': 0.1, 'erev': erev, 'v': -60., 'delay': 0})
    t = st['t']
    p = F.fitpars
    target = st['i']
    # print(F.fitpars)
    pars = F.fit(t, target, F.fitpars)
    fitted = F.exp2syn(st['t'], pars['tau1'], pars['tau2'], pars['weight'], pars['erev'], pars['v'], pars['delay'])
    return(pars, fitted)
    
def plot_all(st, sts=None, fitflag = False):
    P = PH.Plotter((3, 5), figsize=(11, 6))
    ax = P.axarr.ravel()
    fits = {}
    for i, stk in enumerate(st.keys()):
        data = st[stk]
        idat = np.array(data['i'])
        ax[i].plot(data['t'], np.mean(idat, axis=0), 'k-', linewidth=1.5)
        print(idat.shape)
        for j in range(idat.shape[0]):
            ax[i].plot(data['t'], idat[j], 'k-', linewidth=0.5)
        if fitflag:
            fitp, fit = fit_one(data, stk)
            fits[stk] = fitp
            ax[i].plot(data['t'], fit, 'r--', linewidth=0.5)
        ax[i].set_title('%s : %s' % (stk[0], stk[1]), fontsize=7)
    if sts is not None:
        for i, stks in enumerate(sts.keys()):
            datas = sts[stks]
            ax[i].plot(datas['t'], np.mean(datas['i'], axis=0), 'g-', linewidth=3, alpha=0.5)

        # print('fits: ', fits)
    if fitflag:
        stkeys = list(st.keys())
        postcells = ['bushy', 'tstellate', 'dstellate', 'octopus', 'pyramidal', 'tuberculoventral', 'cartwheel']
        fmts = {'weight': '{0:<18.6f}', 'tau1': '{0:<18.3f}', 'tau2': '{0:<18.3f}', 'delay':'{0:<18.3f}', 'erev': '{0:<18.1f}'}
        pars = list(fmts.keys())
        for pre in ['sgc', 'dstellate', 'tuberculoventral']:
            firstline = False
            vrow = dict.fromkeys(pars, False)
            for v in pars:
                for post in postcells:
                    if (pre, post) not in stkeys:
                        print('{0:<18s}'.format('0'), end='')
                        continue
                    if not firstline:
                        print('\n%s' % pre.upper())
                        print('{0:<18s}'.format(' '), end='')
                        for i in range(len(postcells)):
                            print('{0:<18s}'.format(postcells[i]), end='')
                        print()
                        firstline = True
                    if firstline:
                        if not vrow[v]:
                            print('{0:<18s}'.format(v), end='')
                            vrow[v] = True
                    # print(fits[(pre, post)].keys())
                    print(fmts[v].format(fits[(pre, post)][v].value), end='')
                print()

    mpl.show()
    return(st, fits)

def readpkl(mode, mode2=None):
    with(open('%s.pkl' % mode, 'rb')) as fh:
        st = pickle.load(fh)
    st2 = None
    if mode2 is not None:
        with(open('%s.pkl' % mode2, 'rb')) as fh:
            st2 = pickle.load(fh)
    # print('st2: ', st2)
    plot_all(st, st2)
     
def runall(mode='multisite'):
    st = {}
    pkst = {}
    for pre in convergence.keys():
        for post in convergence[pre]:
            # if pre == 'sgc' and post in ['bushy', 'tstellate']:
            sti = runtest([pre, post], mode)
            st[(pre, post)] = sti
            # remove neuron objects
            pkst[(pre, post)] = {'t': sti['t'], 'i': sti.isoma, 'v': sti['v_soma'], 'pre': sti['v_pre']}
    
    with(open('%s.pkl' % mode, 'wb')) as fh:
        pickle.dump(pkst, fh)
    plot_all(pkst)
    

        
def main():
    #testexp()
    #st = runtest()
    mode = 'multisites'
    #st = runall(mode)
    readpkl(mode, 'simple')
    
    if sys.flags.interactive == 0:
        pg.QtGui.QApplication.exec_()
    
if __name__ == '__main__':
    main()  # prevents code from running during sphinx doc generation.
