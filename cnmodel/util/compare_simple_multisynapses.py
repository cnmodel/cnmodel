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

Requires Python 3.6
"""
import sys
import argparse
from pathlib import Path
import numpy as np
# import pyqtgraph as pg
import matplotlib.pyplot as mpl
import cnmodel.util.PlotHelpers as PH
from cnmodel.protocols import SynapseTest
from cnmodel import cells
from cnmodel.synapses import Synapse
import pickle
import lmfit


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
#                                  ('tau2', initpars['tau2'], True, 0.1, 50., None),
                                  ('tauratio', initpars['tauratio'], True, 1.0001, 100., None),
                                  ('weight', initpars['weight'], True, 1e-6, 1, None),
                                  ('erev', initpars['erev'], False), # do not adjust!
                                  ('v', initpars['v'], False),
                                  ('delay', initpars['delay'], True, 0., 5., None))
            self.func = self.exp2syn_err
        else:
            raise ValueError

    def fit(self, x, y, p, verbose=False):
        """
        Perform the curve fit against the specified function
        
        Parameters
        ----------
        x : time base for waveform (np array or list)
        y : waveform (1-d np array or list)
        p : parameters in lmfit Parameters structure
        verbose : boolean (default: False)
            If true, print the parameters in a nice format
        """
        
        kws={'maxfev': 5000}
        # print('p: ', p)
        self.mim = lmfit.minimize(self.func, p, method='least_squares', args=(x, y)) #, kws=kws)
        if verbose:
            lmfit.printfuncs.report_fit(self.mim.params)
        fitpars = self.mim.params
        return fitpars

    #@staticmethod
    def exp2syn(self, x, tau1, tauratio, weight, erev, v, delay):
        """
        Compute the exp2syn waveform current as it is done in Neuron
        
        Note on units:
        The units are assumed to be "self consistent"
        Thus, if the time base is in msec, tau1, tau2 and the delay are
            in msec. erev and v should be in matching units (e.g., mV)
        Note thtat the weight is not equal to the conducdtance G, because
        of the scaling of the waveform
        
        Parameters
        ----------
        x : time base
        tau1 : rising tau for waveform
        tau2 : falling tau for waveform
        weight : amplitude of the waveform (conductance)
        erev : reversal potential for ionic species (used to compute i)
        v : holding voltage at which waveform is computed
        delay : delay to start of function
        
        Returns
        -------
        i, the calucated current trace for these parameters.
        """
        # we handle the requirement that tau2 > tau1 by setting the 
        # expression in lmfit, and using tauratio rather than a direct
        # tau2.
        # if tau1/tau2 > 1.0:  # make sure tau1 is less than tau2
        #     tau1 = 0.999*tau2
        tau2 = tauratio*tau1
        tp = (tau1*tau2)/(tau2 - tau1) * np.log(tau2/tau1)
        factor = -np.exp(-tp/tau1) + np.exp(-tp/tau2)
        factor = 1.0/factor
        G = weight * factor * (np.exp(-(x-delay)/tau2) - np.exp(-(x-delay)/tau1))
        G[x-delay < 0] = 0.
        i = G * (v - erev)  #      i(nanoamps), g(micromhos);

        return i
    
    def exp2syn_err(self, p, x, y):
        return np.fabs(y-self.exp2syn(x, **dict([(k,p.value) for k,p in p.items()])))

    def factor(self, tau1, tauratio, weight):
        """
        calculate tau-scaled weight
        """
        tau2 = tau1*tauratio
        tp = (tau1*tau2)/(tau2 - tau1) * np.log(tau2/tau1)
        factor = -np.exp(-tp/tau1) + np.exp(-tp/tau2)
        factor = 1.0/factor
        G = weight * factor
        return(G)
    
def testexp():
    """
    Test the exp2syn fitting function
    """
    pars = {'tau1': 0.1, 'tauratio': 2.0, 'weight': 0.1, 'erev': -70., 'v': -65., 'delay': 1}
    F = Exp2SynFitting(initpars={'tau1': 0.2,  'tauratio': 2.0, 'weight': 0.1, 'erev': -70., 'v': -65., 'delay': 0})
    t = np.arange(0, 10., 0.01)
    p = F.fitpars
    target = F.exp2syn(t, pars['tau1'], pars['tauratio'], pars['weight'], pars['erev'], pars['v'], pars['delay'])
    # print(F.fitpars)
    pars_fit = F.fit(t, target, F.fitpars)
    print('\nTest fit result: ')
    lmfit.printfuncs.report_fit(F.mim.params)
    print('( tau2 = ', pars_fit['tau1'].value*pars_fit['tauratio'].value, ')')
    print('target parameters: ', pars)
    print('\n')


def compute_psc(synapsetype='multisite', celltypes=['sgc', 'tstellate']):
    """
    Compute the PSC between the specified two cell tpes
    The type of PSC is set by synspase type and must be 'multisite' or 'simple'
    """
    assert(synapsetype in ['multisite', 'simple'])
    
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
    
    print(f"Computing psc for connection {celltypes[0]:s} -> {celltypes[1]:s}")
    nTerminals = convergence.get(celltypes[0], {}).get(celltypes[1], None)
    if nTerminals is None:
        nTerminals = 1
        print(f"Warning: Unknown convergence for {celltypes[0]:s} -> {celltypes[1]:s}, ASSUMING {nTerminals:d} terminals")

    if celltypes == ['sgc', 'bushy']:
        niter = 200
    else:
        niter = 500
    assert(synapsetype in ['simple', 'multisite'])
    st = SynapseTest()
    dt = 0.010
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
    # st.show_result()  # pyqtgraph plotting - 
    # st.plots['VPre'].setYRange(-70., 10.)
    # st.plots['EPSC'].setYRange(-2.0, 0.5)
    # st.plots['latency2080'].setYRange(0., 1.0)
    # st.plots['halfwidth'].setYRange(0., 1.0)
    # st.plots['RT'].setYRange(0., 0.2)
    # st.plots['latency'].setYRange(0., 1.0)
    # st.plots['latency_distribution'].setYRange(0., 1.0)
    return st  # need to keep st alive in memory

def fit_one(st, stk):
    """
    Fit one trace to the exp2syn function
    
    Parameters
    ----------
    
    st : dict (no default)
        A dictionary containing (at least) the following keys:
        't' : the time base for the trace to be fit
        'i' : a list of arrays or a 2d numpy array of the data to be fit
    
    Returns
    -------
    pars : The fitting parameters 
        The fitted parameters are returned as an lmfit Parameters object, so access
        individual parameters a pars['parname'].value
    
    fitted : The fitted waveform, of the same length as the source data waveform
    """
    

    if stk[0] in ['sgc', 'tstellate', 'granule']:
        erev = 0.  # set erev according to the source cell (excitatory: 0, inhibitory: -80)
    else:
        erev = -70.  # value used in gly_psd
    print('\nstk, erev: ', stk, erev)
    F = Exp2SynFitting(initpars={'tau1': 1.0, 'tauratio': 4.0, 'weight': 0.001, 'erev': erev, 'v': -65., 'delay': 0})
    t = st['t']
    p = F.fitpars
    target = np.mean(np.array(st['i']), axis=0)
    # print(F.fitpars)
    pars = F.fit(t, target, F.fitpars)
    gw = F.factor(pars['tau1'].value, pars['tauratio'].value, pars['weight'].value)
    pars.add('GWeight', value=gw, vary=False)
    fitted = F.exp2syn(st['t'], pars['tau1'], pars['tauratio'], pars['weight'], pars['erev'], pars['v'], pars['delay'])
    lmfit.printfuncs.report_fit(F.mim.params)
    return(pars, fitted)


def fit_all():
    """
    Fit exp2syn against the traces in stm
    This fits all pre-post cell pairs, and returns the fit
    
    """
    stm = read_pickle('multisite.pkl')
    fits = {}
    fitted = {}
    # print('stm keys: ', stm.keys())
    # exit()
    for i, stk in enumerate(stm.keys()): # for each pre-post cell pair
        fitp, fit = fit_one(stm[stk], stk)
        fits[stk] = fitp

        fitted[stk] = {'t': stm[stk]['t'], 'i': [fit], 'pars': fitp}
    with(open('simple.pkl', 'wb')) as fh:
        pickle.dump(fitted, fh) 
    return(fits)       
        
        
def plot_all(stm, sts=None):
    P = PH.Plotter((3, 5), figsize=(11, 6))
    ax = P.axarr.ravel()
    keypairorder = []
    for i, stk in enumerate(stm.keys()):
        keypairorder.append(stk)
        data = stm[stk]
        idat = np.array(data['i'])
        ax[i].plot(np.array(data['t']), np.mean(idat, axis=0), 'c-', linewidth=1.5)
        sd = np.std(idat, axis=0)
        ax[i].plot(np.array(data['t']), np.mean(idat, axis=0)+sd, 'c--', linewidth=0.5)        # for j in range(idat.shape[0]):
        ax[i].plot(np.array(data['t']), np.mean(idat, axis=0)-sd, 'c--', linewidth=0.5)        # for j in range(idat.shape[0]):
        rel = 0
        for j in range(idat.shape[0]):
            if np.min(idat[j,:]) < -1e-2 or np.max(idat[j,:]) > 1e-2:
                rel += 1
        print(f'{str(stk):s} {rel:d} of {idat.shape[0]:d} {str(idat.shape):s}')
            
        #     ax[i].plot(data['t'], idat[j], 'k-', linewidth=0.5, alpha=0.25)
        ax[i].set_title('%s : %s' % (stk[0], stk[1]), fontsize=7)
    
    if sts is not None:  # plot the matching exp2syn on this
        for i, stks in enumerate(keypairorder):
            datas = sts[stks]
            isdat = np.array(datas['i'])
            ax[i].plot(datas['t'], np.mean(isdat, axis=0), 'm-', linewidth=3, alpha=0.5)
    mpl.show()


def print_exp2syn_fits(st):
    stkeys = list(st.keys())
    postcells = ['bushy', 'tstellate', 'dstellate', 'octopus', 'pyramidal', 'tuberculoventral', 'cartwheel']
    fmts = {'weight': '{0:<18.6f}', 'GWeight': '{0:<18.6f}', 'tau1': '{0:<18.3f}', 'tau2': '{0:<18.3f}', 'delay':'{0:<18.3f}', 'erev': '{0:<18.1f}'}
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
                        if v == 'tauratio':
                            print('{0:<18s}'.format('tau2'), end='')
                        else:
                            print('{0:<18s}'.format(v), end='')
                        vrow[v] = True
                # print(fits[(pre, post)].keys())
                fits = st[(pre, post)]['pars']
                # print(v, fits)
                if v in ['tauratio', 'tau2']:
                    print(fmts[v].format(fits['tau1'].value*fits['tauratio'].value), end='')
                else:
                    print(fmts[v].format(fits[v].value), end='')
            print()


def read_pickle(mode):
    """
    Read the pickled file generated by runall
    
    Parameters
    ----------
    mode : str
        either 'simple' or 'multisite'
    
    Returns:
        the resulting data, which is a dictionary
    """
    with(open(Path(mode).with_suffix('.pkl'), 'rb')) as fh:
        st = pickle.load(fh)
    return(st)


def run_all():
    """
    Run all of the multisite synapse calculations for each cell pair
    Save the results in the multisite.pkl file
    """
    st = {}
    pkst = {}
    for pre in convergence.keys():
        for post in convergence[pre]:
            # if pre == 'sgc' and post in ['bushy', 'tstellate']:
            sti = compute_psc(synapsetype='multisite', celltypes=[pre, post])
            st[(pre, post)] = sti
            # remove neuron objects before pickling
            pkst[(pre, post)] = {'t': sti['t'], 'i': sti.isoma, 'v': sti['v_soma'], 'pre': sti['v_pre']}
    
    with(open(Path('multisite').with_suffix('.pkl'), 'wb')) as fh:
        pickle.dump(pkst, fh)
    plot_all(pkst)


def main():
    parser = argparse.ArgumentParser(description='Compare simple and multisite synapses')
    parser.add_argument('-m', '--mode', type=str, dest='mode', default='None',
                        choices = ['simple', 'multisite', 'both'],
                        help='Select mode [simple, multisite, compare]')
    parser.add_argument('-r', '--run', action='store_true', dest='run',
                        help='Run multisite models between all cell types')
    parser.add_argument('-f', '--fit', action='store_true', dest='fit',
                        help='Fit exp2syn waveforms to multisite data')
    parser.add_argument('-t', '--test', action='store_true', dest='test',
                        help='Run the test on exp2syn fitting')
    parser.add_argument('-p', '--plot', action='store_true', dest='plot',
                        help='Plot the current comparison between simple and multisite')
    parser.add_argument('-l', '--list', action='store_true', dest='list',
                        help='List the simple fit parameters')

    args = parser.parse_args()
    
    if args.test:
        testexp()
        exit()

    if args.run:
        run_all()

    if args.fit:
        fit_all() # self contained - always fits exp2syn against the current multistie data
        ds = read_pickle('simple')
        dm = read_pickle('multisite')
        plot_all(dm, ds)
        exit()
        
    if args.plot:
        if args.mode in ['simple', 'multisite']:
            d = read_pickle(args.mode)
            plot_all(d)
            exit()
        elif args.mode in ['both']:
            ds = read_pickle('simple')
            dm = read_pickle('multisite')
            plot_all(dm, ds)
            exit()
        else:
            print(f"Mode {args.mode:s} is not valid")
            exit()
            

    if args.list:
        ds = read_pickle('simple')
        print_exp2syn_fits(ds)
    
    # if sys.flags.interactive == 0:
    #     pg.QtGui.QApplication.exec_()
    
if __name__ == '__main__':
    main()  
