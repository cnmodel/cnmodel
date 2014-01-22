#!/usr/bin/python
#
# Synapse definitions for models.
#
# This file includes a number of different synapse definitions and default conductances
# for point models. Most are models from the lab for neurons of the cochlear nucleus.
# the synaptic receptor models are gleaned from the literature and sometimes fitted to the 
# cochlear nucleus data.
#
# Paul B. Manis, Ph.D. 2009 (August - November 2009)
#

import os, os.path
import gc
from optparse import OptionParser
from itertools import chain
from neuron import h
from neuron import *
from math import sqrt
from math import pi as PI
import numpy
import scipy
import scipy.interpolate
from scipy.stats import gamma
from numpy.random import lognormal
from nrnlibrary.pynrnutilities import *
import nrnlibrary.Cells as Cells
import matplotlib.pylab as mpl
import pdb

mpl.rcParams['interactive'] = False

runQuiet = True
veryFirst = 1
# utility class to create parameter lists... 
# create like: p = Params(abc=2.0, defg = 3.0, lunch='sandwich')
# reference like p.abc, p.defg, etc.
class Params(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


def template_SGAxon(debug=False, message=None):
    axon = Cells.hh(message='')
    return (axon)

# define the "Calyx" template.
# Each axon ends in (and thus includes) a calyx; therefore we make a list of axons
# based on the template.
# incoming argument: which axon (index) to address SGs
# Technically, this isn't the axon, but the terminal....

def template_Calyx_Billup(debug=False, nzones=1, message=None):
    calyx = Cells.hh(message='  >> creating calyx_Billup')
    coh = []
    calyx.push()
    for k in range(0, nzones):
        coh.append(h.COH2(0.5, calyx))
    h.pop_section()
    return (calyx, coh)


def template_multisite(debug=False, nzones=1, celltype='bushy', message=None,
                       type='lognormal', Identifier=0, stochasticPars=None, calciumPars=None):
    if stochasticPars is None:
        raise TypeError
    global veryFirst
    mu = u'\u03bc'
    sigma = u'\u03c3'
    term = Cells.hh(
        message='  >> creating terminal with %d release zones using lognormal release latencies (coh4)' % nzones)
    terminal = term[0]
    terminal.push()
    cleft = []
    if calciumPars is not None:
        terminal.insert('cap') # insert calcium channel density
        terminal().cap.pcabar = calciumPars.Ca_gbar
        terminal.insert('cad')
    relsite = h.COH4(0.5, sec=terminal)
    relsite.nZones = nzones
    relsite.rseed = 2 # int(numpy.random.random_integers(1,1024))
    relsite.latency = stochasticPars.latency
    relsite.latstd = stochasticPars.LN_std
    if debug is True:
        relsite.debug = 1
    relsite.Identifier = Identifier
    # if type == 'gamma':
    #     gd = gamma.rvs(2, size=10000)/2.0 # get a sample of 10000 events with a gamma dist of 2, set to mean of 1.0
    #     if relsite.latstd > 0.0:
    #         gds = relsite.latency+std*(gd-1.0)/gd.std() # scale standard deviation
    #     else:
    #         gds = relsite.latency*numpy.ones((10000,1))
    # if type == 'lognormal':
    #     if std > 0.0:
    #         gds = lognormal(mean=0, sigma=relsite.latstd, size=10000)
    #     else:
    #         gds = numpy.zeros((10000, 1))
    # use the variable latency mode of COH4. And, it is lognormal no matter what.
    # the parameters are defined in COH4.mod as follows
    # 	 Time course of latency shift in release during repetitive stimulation
    # 	Lat_Flag = 0 (1) : 0 means fixed latency, 1 means lognormal distribution
    # 	Lat_t0 = 0.0 (ms) : minimum time since simulation start before changes in latency are calculated
    # 	Lat_A0 = 0.0 (ms) : size of latency shift from t0 to infinity
    # 	Lat_tau = 100.0 (ms) : rate of change of latency shift (from fit of a+b(1-exp(-t/tau)))
    # 	: Statistical control of log-normal release shape over time during repetive stimulation
    # 	LN_Flag = 0 (1) : 0 means fixed values for all time
    # 	LN_t0 = 0.0 (ms) : : minimum time since simulation start before changes in distribution are calculated
    # 	LN_A0 = 0.0 (ms) : size of change in sigma from t0 to infinity
    # 	LN_tau = 100.0 (ms) : rate of change of sigma over time (from fit of a+b*(1-exp(-t/tau)))

    relsite.LN_Flag = stochasticPars.LN_Flag # enable use of lognormal release latency
    relsite.LN_t0 = stochasticPars.LN_t0
    relsite.LN_A0 = stochasticPars.LN_A0
    relsite.LN_tau = stochasticPars.LN_tau
    relsite.Lat_Flag = stochasticPars.Lat_Flag
    relsite.Lat_t0 = stochasticPars.Lat_t0
    relsite.Lat_A0 = stochasticPars.Lat_A0
    relsite.Lat_tau = stochasticPars.Lat_tau
    if veryFirst == 1 and debug is True:
        veryFirst = 0
        #mpl.figure(2)
    if celltype == 'bushy':
        relsite.TDur = 0.10
        relsite.TAmp = 0.770
    else:
        relsite.TDur = 0.25
        relsite.TAmp = 1.56625
    h.pop_section()
    if not runQuiet:
        print "adding %d clefts to terminal " % nzones, terminal
    for k in range(0, nzones):
        terminal.push()
        cl = h.cleftXmtr(0.5, sec=terminal)
        h.pop_section()
        cleft.append(cl) # cleft
    return (terminal, relsite, cleft)

################################################################################
# define the "PSD" template.
# Each PSD has receptors for each active zone, matched to presynaptic
# terminals.
# Template requires a segment to insert the receptors into as first argument
# also requires a number of receptor sites to insert as second argument
# e.g., psd = new PSD(section, nReceptors)
################################################################################

def template_iGluR_PSD(sec, nReceptors=1, debug=False, cellname=None, message=None, NMDARatio=1):
    psd = []
    psdn = []
    sec.push()
    AN_Po_Ratio = 23.2917 # ratio of open probabilities for AMPA and NMDAR's at peak currents
    for k in range(0, nReceptors):
        psd.append(h.AMPATRUSSELL(0.5, sec)) # raman/trussell AMPA with rectification
        psdn.append(h.NMDA_Kampa(0.5, sec)) # Kampa state model NMDA receptors

        if cellname == 'bushy':
            psd[-1].Ro1 = 107.85
            psd[-1].Ro2 = 0.6193
            psd[-1].Rc1 = 3.678
            psd[-1].Rc2 = 0.3212
            psdn[-1].gNAR = 0.036 * AN_Po_Ratio * NMDARatio # 0.36*AN_Po_Ratio*NMDARatio
            if k == 0:
                print "Bushy NMDAR set to %8.2f" % psdn[-1].gNAR
        if cellname == 'stellate':
            psd[-1].Ro1 = 39.25
            psd[-1].Ro2 = 4.40
            psd[-1].Rc1 = 0.667
            psd[-1].Rc2 = 0.237
            psd[-1].PA = 0.1
            psdn[-1].gNAR = 1 * AN_Po_Ratio * NMDARatio

    h.pop_section()
    par = {'Ro1': ('r', psd[0].Ro1),
           'Ro2': ('r', psd[0].Ro2),
           'Rc1': ('r', psd[0].Rc1),
           'Rc2': ('r', psd[0].Rc2), }
    parn = {'gNAR': ('4', psdn[0].gNAR), }
    return (psd, psdn, par, parn)


# the following templates are a bit more complicated.
# The parameter names as defined in the model are returned
# but also those that are involved in the forward binding reactions are
# listed separately - this is to allow us to run curve fits that adjust
# only a subset of the parameters a  time - e.g., the rising phase, then
# the falling phase with the rising phase fixed.
# the dictionary selection is made by selectpars in glycine_fit.py.
#

def template_Gly_PSD_exp(sec, debug=False, nReceptors=2, cellname=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLY2(0.5, sec))
    h.pop_section()
    par = ['alpha', 'beta']
    p = []
    for n in par:
        p.append(eval('psd[0].' + n)) # get default values from the mod file
    return (psd, par, p)


def template_Gly_PSD_State_Glya5(sec, debug=False, nReceptors=2, psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLYa5(0.5, sec))
    h.pop_section()
    par = {'kf1': ('r', psd[0].kf1), # retreive values in the MOD file
           'kf2': ('r', psd[0].kf2),
           'kb1': ('r', psd[0].kb1),
           'kb2': ('r', psd[0].kb2),
           'a1': ('f', psd[0].a1),
           'b1': ('f', psd[0].b1),
           'a2': ('f', psd[0].a2),
           'b2': ('f', psd[0].b2)}
    return (psd, par)


def template_Gly_PSD_State_Gly6S(sec, debug=False, nReceptors=2, psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.Gly6S(0.5, sec)) # simple using Trussell model 6 states with desens
        if not runQuiet:
            print 'Gly6S psdtype: ', psdtype
        if psdtype == 'glyslow': # fit on 8 March 2010, error =  0.164, max open: 0.155
            psd[-1].Rd = 1.177999
            psd[-1].Rr = 0.000005
            psd[-1].Rb = 0.009403
            psd[-1].Ru2 = 0.000086
            psd[-1].Ro1 = 0.187858
            psd[-1].Ro2 = 1.064426
            psd[-1].Ru1 = 0.028696
            psd[-1].Rc1 = 0.103625
            psd[-1].Rc2 = 1.730578
    h.pop_section()
    par = {'Rb': ('n', psd[0].Rb), # retrive values in the MOD file
           'Ru1': ('r', psd[0].Ru1),
           'Ru2': ('r', psd[0].Ru2),
           'Rd': ('f', psd[0].Rd),
           'Rr': ('f', psd[0].Rr),
           'Ro1': ('f', psd[0].Ro1),
           'Ro2': ('r', psd[0].Ro2),
           'Rc1': ('f', psd[0].Rc1),
           'Rc2': ('r', psd[0].Rc2)}
    return (psd, par)


def template_Gly_PSD_State_PL(sec, debug=False, nReceptors=2, cellname=None,
                              psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLYaPL(0.5, sec)) # simple dextesche glycine receptors
        if not runQuiet:
            print 'PL psdtype: ', psdtype
        if psdtype == 'glyslow':
            psd[-1].a1 = 0.000451
            psd[-1].a2 = 0.220
            psd[-1].b1 = 13.27
            psd[-1].b2 = 6.845
            psd[-1].kon = 0.00555
            psd[-1].koff = 2.256
            psd[-1].r = 1.060
            psd[-1].d = 55.03
        if psdtype == 'glyfast': # fit from 3/5/2010. error =  0.174 maxopen = 0.0385
            psd[-1].a1 = 1.000476
            psd[-1].a2 = 0.137903
            psd[-1].b1 = 1.700306
            psd[-1].koff = 13.143132
            psd[-1].kon = 0.038634
            psd[-1].r = 0.842504
            psd[-1].b2 = 8.051435
            psd[-1].d = 12.821820
    h.pop_section()
    par = {'kon': ('r', psd[0].kon), # retrive values in the MOD file
           'koff': ('r', psd[0].koff),
           'a1': ('r', psd[0].a1),
           'b1': ('r', psd[0].b1),
           'a2': ('f', psd[0].a2),
           'b2': ('f', psd[0].b2),
           'r': ('f', psd[0].r),
           'd': ('f', psd[0].d)}
    return (psd, par)


def template_Gly_PSD_State_GC(sec, debug=False, nReceptors=2, psdtype=None, message=None):
    psd = []
    sec.push()
    for k in range(0, nReceptors):
        psd.append(h.GLYaGC(0.5, sec)) # simple dextesche glycine receptors
        if psdtype == 'glyslow':
            psd[-1].k1 = 12.81    # (/uM /ms)	: binding
            psd[-1].km1 = 0.0087#	(/ms)	: unbinding
            psd[-1].a1 = 0.0195#	(/ms)	: opening
            psd[-1].b1 = 1.138    #(/ms)	: closing
            psd[-1].r1 = 6.13    #(/ms)	: desense 1
            psd[-1].d1 = 0.000462 # (/ms)	: return from d1
            psd[-1].r2 = 0.731# (/ms)     : return from deep state
            psd[-1].d2 = 1.65 #(/ms)     : going to deep state
            psd[-1].r3 = 3.83 #(/ms)     : return from deep state
            psd[-1].d3 = 1.806 #(/ms)     : going to deep state
            psd[-1].rd = 1.04 #(/ms)
            psd[-1].dd = 1.004 # (/ms)
    h.pop_section()
    par = {'k1': ('r', psd[0].k1), # retrive values in the MOD file
           'km1': ('r', psd[0].km1),
           'a1': ('r', psd[0].a1),
           'b1': ('r', psd[0].b1),
           'r1': ('f', psd[0].r1),
           'd1': ('f', psd[0].d1),
           'r2': ('f', psd[0].r2),
           'd2': ('f', psd[0].d2),
           'r3': ('f', psd[0].r3),
           'd3': ('f', psd[0].d3),
           'rd': ('f', psd[0].rd),
           'dd': ('f', psd[0].dd)}
    return (psd, par)


################################################################################
# The following routines set the synapse dynamics, based on measurements and fit
# to the Dittman-Regehr model.
################################################################################

def setDF(coh, celltype, synapsetype, select=None):
    """ set the parameters for the calyx release model ...
        These paramteres were obtained from an optimized fit of the Dittman-Regehr
        model to stimulus and recovery data for the synapses at 100, 200 and 300 Hz,
        for times out to about 0.5 - 1.0 second. Data from Ruili Xie and Yong Wang.
        Fitting by Paul Manis
    """
    if celltype == 'bushy':
        if synapsetype == 'epsc':
            coh = bushy_epsc(coh)
        if synapsetype == 'ipsc':
            if select is None:
                coh = bushy_ipsc_average(coh)
            else:
                coh = bushy_ipsc_single(coh, select=select)
    if celltype == 'stellate':
        if synapsetype == 'epsc':
            coh = stellate_epsc(coh)
        if synapsetype == 'ipsc':
            coh = stellate_ipsc(coh)
    return (coh)


def bushy_epsc(synapse):
    """ data is average of 3 cells studied with recovery curves and individually fit """
    synapse.F = 0.29366
    synapse.k0 = 0.52313 / 1000.0
    synapse.kmax = 19.33805 / 1000.0
    synapse.taud = 15.16
    synapse.kd = 0.11283
    synapse.taus = 17912.2
    synapse.ks = 11.531
    synapse.kf = 17.78
    synapse.tauf = 9.75
    synapse.dD = 0.57771
    synapse.dF = 0.60364
    synapse.glu = 2.12827
    return (synapse)


def stellate_epsc(synapse):
    """ data is average of 3 cells studied with recovery curves and individually fit """
    synapse.F = 0.43435
    synapse.k0 = 0.06717 / 1000.0
    synapse.kmax = 52.82713 / 1000.0
    synapse.taud = 3.98
    synapse.kd = 0.08209
    synapse.taus = 16917.120
    synapse.ks = 14.24460
    synapse.kf = 18.16292
    synapse.tauf = 11.38
    synapse.dD = 2.46535
    synapse.dF = 1.44543
    synapse.glu = 5.86564
    return (synapse)


def stellate_ipsc(synapse):
    """ data is average of 3 cells studied with recovery curves and individually fit, 100 Hz """
    synapse.F = 0.23047
    synapse.k0 = 1.23636 / 1000.0
    synapse.kmax = 45.34474 / 1000.0
    synapse.taud = 98.09
    synapse.kd = 0.01183
    synapse.taus = 17614.50
    synapse.ks = 17.88618
    synapse.kf = 19.11424
    synapse.tauf = 32.28
    synapse.dD = 2.52072
    synapse.dF = 2.33317
    synapse.glu = 3.06948
    return (synapse)


def bushy_ipsc_average(synapse):
    """average of 16 Bushy cells. Done differently than other averages.
    The individual fits were compiled, and an average computed for just the 100 Hz data
    across the individual fits. This average was then fit to the function of Dittman and Regeher
    (also in Xu-Friedman's papers). 
    The individual cells show a great deal of variability, from straight depression, to 
    depression/facilitaiton mixed, to facilation alone. This set of parameters generates
    a weak facilitaiton followed by depression back to baseline.
    """
    print "USING average kinetics for Bushy IPSCs"

    # average of 16cells for 100 Hz (to model); no recovery.
    synapse.F = 0.18521
    synapse.k0 = 2.29700
    synapse.kmax = 27.6667
    synapse.taud = 0.12366
    synapse.kd = 0.12272
    synapse.taus = 9.59624
    synapse.ks = 8.854469
    synapse.kf = 5.70771
    synapse.tauf = 0.37752
    synapse.dD = 4.00335
    synapse.dF = 0.72605
    synapse.glu = 5.61985

    # estimates for 400 Hz
    # synapse.F = 0.09
    # synapse.k0 = 1.2;
    # synapse.kmax = 30.;
    # synapse.taud = 0.01
    # synapse.kd = 0.75
    # synapse.taus = 0.015
    # synapse.ks = 1000.0
    # synapse.kf = 5.0
    # synapse.tauf = 0.3
    # synapse.dD = 1.0
    # synapse.dF = 0.025
    # synapse.glu = 4

    # synapse.F = 0.085426
    # synapse.k0 = 1.199372
    # synapse.kmax = 24.204277
    # synapse.taud = 0.300000
    # synapse.kd = 1.965292
    # synapse.taus = 2.596443
    # synapse.ks = 0.056385
    # synapse.kf = 0.721157
    # synapse.tauf = 0.034560
    # synapse.dD = 0.733980
    # synapse.dF = 0.025101
    # synapse.glu = 3.877192

    # average of 8 cells, all at 100 Hz with no recovery
    # synapse.F =    0.2450
    # synapse.k0 =   1.6206/1000.0
    # synapse.kmax = 26.0607/1000.0
    # synapse.taud = 0.0798
    # synapse.kd =   0.9679
    # synapse.taus = 9.3612
    # synapse.ks =   14.3474
    # synapse.kf =    4.2168
    # synapse.tauf =  0.1250
    # synapse.dD =    4.2715
    # synapse.dF =    0.6322
    # synapse.glu =   8.6160

    # average of 5 cells, mostly 100 Hz, but some 50, 200 and 400
    # synapse.F =    0.15573
    # synapse.k0 =   2.32272/1000.
    # synapse.kmax = 28.98878/1000.
    # synapse.taud = 0.16284
    # synapse.kd =   2.52092
    # synapse.taus = 17.97092
    # synapse.ks =   19.63906
    # synapse.kf =   7.44154
    # synapse.tauf = 0.10193
    # synapse.dD =   2.36659
    # synapse.dF =   0.38516
    # synapse.glu =  8.82600

    #original average - probably skewed.
    # synapse.F = 0.23382
    # synapse.k0 = 0.67554/1000.0
    # synapse.kmax = 52.93832/1000.0
    # synapse.taud = 8.195
    # synapse.kd = 0.28734
    # synapse.taus = 17.500
    # synapse.ks = 4.57098
    # synapse.kf = 16.21564
    # synapse.tauf = 123.36
    # synapse.dD = 2.21580
    # synapse.dF = 1.17146
    # synapse.glu = 1.90428
    return (synapse)


def bushy_ipsc_single(synapse, select=None):
    """ data is from 31aug08b (single cell, clean dataset)"""
    print "Using bushy ipsc"

    if select is None or select > 4 or select <= 0:
        return (bushy_ipsc_average(synapse))

    if select is 1: # 30aug08f
        print "using 30aug08f ipsc"
        synapse.F = 0.221818
        synapse.k0 = 0.003636364
        synapse.kmax = 0.077562107
        synapse.taud = 0.300000
        synapse.kd = 1.112554
        synapse.taus = 3.500000
        synapse.ks = 0.600000
        synapse.kf = 3.730452
        synapse.tauf = 0.592129
        synapse.dD = 0.755537
        synapse.dF = 2.931578
        synapse.glu = 1.000000

    if select is 2: #30aug08h
        print "using 30aug08H ipsc"
        synapse.F = 0.239404
        synapse.k0 = 3.636364 / 1000.
        synapse.kmax = 16.725479 / 1000.
        synapse.taud = 0.137832
        synapse.kd = 0.000900
        synapse.taus = 3.500000
        synapse.ks = 0.600000
        synapse.kf = 4.311995
        synapse.tauf = 0.014630
        synapse.dD = 3.326148
        synapse.dF = 0.725512
        synapse.glu = 1.000000

    if select is 3:
        print "using IPSC#3 "
        synapse.F = 0.29594
        synapse.k0 = 0.44388 / 1000.0
        synapse.kmax = 15.11385 / 1000.0
        synapse.taud = 0.00260
        synapse.kd = 0.00090
        synapse.taus = 11.40577
        synapse.ks = 27.98783
        synapse.kf = 30.00000
        synapse.tauf = 0.29853
        synapse.dD = 3.70000
        synapse.dF = 2.71163
        synapse.glu = 4.97494

    return (synapse)


# make the calyx synapses (or just endings on stellate cells)
def stochastic_synapses(h, axon=None, targetcell=None, nFibers=1, nRZones=1,
                        cellname='bushy', psdtype='ampa', message=None, debug=False,
                        thresh=0.0, gmax=1.0, gvar=0, eRev=0,
                        stochasticPars=None, calciumPars=None,
                        NMDARatio=1.0, select=None, Identifier=0):
    """ This routine generates the synaptic connections from "nFibers" presynaptic
        inputs onto a postsynaptic cell.
        Each connection is a stochastic presynaptic synapse ("presynaptic") with
        depression and facilitation.
        Each synapse can have  multiple releasesites "releasesites" on the
        target cell, as set by "NRZones".
        Each release site releases transmitter using "cleftXmtr"
        Each release site's transmitter is in turn attached to a PSD at each ending ("psd")
        Each psd can have a different conductance centered about the mean of
        gmax, according to a gaussian distribution set by gvar.
        NOTE: stochasticPars must define parameters used by multisite, including:
            .delay is the netcon delay between the presynaptic AP and the start of release events
            .Latency is the latency to the mean release event... this could be confusing.
    """
    if stochasticPars is None:
        raise TypeError
    if axon is None:
        raise TypeError
        exit()
        #if debug is True:
    print "\nTarget cell  = %s, psdtype = %s" % (cellname, psdtype)
    printParams(stochasticPars)
    glyslowPoMax = 0.162297  # thse were measured from the kinetic models in Synapses.py, as peak open P for the glycine receptors
    glyfastPoMax = 0.038475  # also later verified, same numbers...
    if psdtype == 'glyfast':
        gmax /= glyfastPoMax  # normalized to maximum open probability for this receptor
    if psdtype == 'glyslow':
        gmax /= glyslowPoMax  # normalized to max open prob for the slow receptor.
    presynaptic = []
    releasesites = []
    allclefts = []
    allpsd = []
    allpar = []
    netcons = [] # build list of connections from individual release sites to the mother calyx
    if isinstance(axon, tuple):
        axon = axon[0]
    for j in range(0, nFibers):  # for each input fiber to the target cell
    #        print "Stochastic syn: j = %d of nFibers = %d nRZones = %d\n" % (j, nFibers, nRZones)
        # for each fiber, create a presynaptic ending with nzones release sites
        (terminal, relzone, clefts) = template_multisite(nzones=nRZones, celltype=cellname,
                                                         message='   synapse %d' % j,
                                                         stochasticPars=stochasticPars,
                                                         calciumPars=calciumPars,
                                                         Identifier=Identifier, debug=debug)
        # and then make a set of postsynaptic zones on the postsynaptic side
        #        print 'PSDTYPE: ', psdtype
        if psdtype == 'ampa':
            relzone = setDF(relzone, cellname, 'epsc') # set the parameters for release
        elif psdtype == 'glyslow' or psdtype == 'glyfast' or psdtype == 'glyexp' or psdtype == 'glya5' or psdtype == 'glyGC':
            relzone = setDF(relzone, cellname, 'ipsc', select) # set the parameters for release
        if psdtype == 'ampa':
            (psd, psdn, par, parn) = template_iGluR_PSD(targetcell, nReceptors=nRZones,
                                                        cellname=cellname, NMDARatio=NMDARatio)
        elif psdtype == 'glyslow':
            (psd, par) = template_Gly_PSD_State_Gly6S(targetcell, nReceptors=nRZones,
                                                      psdtype=psdtype)
        elif psdtype == 'glyfast':
            (psd, par) = template_Gly_PSD_State_PL(targetcell, nReceptors=nRZones,
                                                   psdtype=psdtype)
        elif psdtype == 'glyGC':
            (psd, par) = template_Gly_PSD_State_GC(targetcell, nReceptors=nRZones,
                                                   psdtype=psdtype)
        elif psdtype == 'glya5':
            (psd, par) = template_Gly_PSD_State_Glya5(targetcell, nReceptors=nRZones,
                                                      psdtype=psdtype)
        elif psdtype == 'glyexp':
            (psd, par) = template_Gly_PSD_exp(targetcell, nReceptors=nRZones,
                                              psdtype=psdtype)
        else:
            print "**PSDTYPE IS NOT RECOGNIZED: [%s]\n" % (psdtype)
            exit()
            # connect the mechanisms on the presynaptic side together
        print 'terminal: ', terminal
        print 'axon: ', axon
        terminal.connect(axon, 1, 0) # 1. connect the axon to the current calyx
        terminal.push()
        netcons.append(h.NetCon(terminal(0.5)._ref_v, relzone, thresh, stochasticPars.delay, 1.0))
        netcons[-1].weight[0] = 1
        netcons[-1].threshold = -30.0
        # ****************
        # delay CV is about 0.3 (Sakaba and Takahashi, 2001), assuming delay is about 0.75
        # This only adjusts the delay to each "terminal" from a single an - so is path length,
        # but is not the individual site release latency. That is handled in COH4
        #        if cellname == 'bushy':
        #            delcv = 0.3*(delay/0.75) # scale by specified delay
        #            newdelay = delay+delcv*numpy.random.standard_normal()
        #            print "delay: %f   newdelay: %f   delcv: %f" % (delay, newdelay, delcv)
        #            netcons[-1].delay = newdelay # assign a delay to EACH zone that is different...
        #*****************

        for k in range(0, nRZones): # 2. connect each release site to the mother axon
            if psdtype == 'ampa': # direct connection, no cleft
                relzone.setpointer(relzone._ref_XMTR[k], 'XMTR', psd[k])
                relzone.setpointer(relzone._ref_XMTR[k], 'XMTR', psdn[k]) # include NMDAR's as well at same release site
            elif psdtype == 'glyslow' or psdtype == 'glyfast' or psdtype == 'glyexp' or psdtype == 'glya5' or psdtype == 'glyGC':
                relzone.setpointer(relzone._ref_XMTR[k], 'pre',
                                   clefts[k]) # connect the cleft to the release of transmitter
                clefts[k].preThresh = 0.1
                if cellname == 'stellate' or cellname == 'test':
                    clefts[k].KV = 531.0 # set cleft transmitter kinetic parameters
                    clefts[k].KU = 4.17
                    clefts[k].XMax = 0.731
                elif cellname == 'bushy':
                    clefts[k].KV = 1e9 # really fast for bushy cells.
                    clefts[k].KU = 4.46
                    clefts[k].XMax = 0.733
                else:
                    print('Error in cell name, dying:   '), cellname
                    exit()
                clefts[k].setpointer(clefts[k]._ref_CXmtr, 'XMTR', psd[k]) #connect transmitter release to the PSD
            else:
                print "PSDTYPE IS NOT RECOGNIZED: [%s]\n" % (psdtype)
                exit()
            v = 1.0 + gvar * numpy.random.standard_normal()
            psd[k].gmax = gmax * v # add a little variability - gvar is CV of amplitudes
            psd[k].Erev = eRev # set the reversal potential
            if psdtype == 'ampa': # also adjust the nmda receptors at the same synapse
                psdn[k].gmax = gmax * v
                psdn[k].Erev = eRev
        h.pop_section()
        presynaptic.append(terminal) # add terminal to the presynaptic list of terminals
        releasesites.append(relzone) # build the release site list
        allpsd.extend(psd) # keep a list of the psds as well
        allpar.extend(par)
        if psdtype == 'ampa':
            allpsd.extend(psdn)
            allpar.extend(parn)
        if psdtype is not 'ampa' and psdtype is not 'nmda': # we don't use "cleft" in the AMPA model
            allclefts.extend(clefts)
        if not runQuiet:
            if message is not None:
                print message
    return (presynaptic, releasesites, allpsd, allclefts, netcons, allpar)


def print_model_parameters(psd, par, message=None):
    if message is not None:
        print "%s" % message
    for i in range(0, len(par)):
        for j in range(0, len(psd[i])):
            print '%s: ' % (psd[i].hname()),
            for k in par[i]:
                v = eval('psd[i][j].' + k)
                print '%8s  = %12f' % (k, v)


def ParameterInit(delay=0.6, latency=0.0, std=0.0):
    """ convenience function to set up a default parameter group for stochastic synapses, in which
    there is no time variance in release or release latency
    Note the following: these parameters are passed to coh4.mod. 
    In that routine, if LN_Flag is 0, then there is no time evolution of the distribution,
    BUT the sigma value is passed (called "std") and is used to jitter the release events.
    If you don't want this, set LN_std to 0.0
    """
    vPars = Params(LN_Flag=0, LN_t0=10.0, LN_A0=0., LN_tau=10.0, LN_std=std,
                   Lat_Flag=0, Lat_t0=10.0, Lat_A0=0., Lat_tau=10.,
                   delay=delay, latency=latency)
    return (vPars)


def printParams(p):
    """
    print the parameter block created in Parameter Init"""
    u = dir(p)
    print "-------- Synapse Parameter Block----------"
    for k in u:
        if k[0] != '_':
            print "%15s = " % (k), eval('p.%s' % k)
    print "-------- ---------------------- ----------"


def make_pulse(stim, pulsetype='square'):
    global j
    delay = int(numpy.floor(stim['delay'] / h.dt))
    ipi = int(numpy.floor((1000.0 / stim['Sfreq']) / h.dt))
    pdur = int(numpy.floor(stim['dur'] / h.dt))
    posttest = int(numpy.floor(stim['PT'] / h.dt))
    maxt = h.dt * (stim['delay'] + (ipi * (stim['NP'] + 3)) + posttest + pdur * 5)
    w = numpy.zeros(floor(maxt / h.dt))
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
                    w[delay + ipi * j + i] += stim['amp'] * (1.0 - exp(-i / (pdur / 3.0))) * exp(
                        -(i - (pdur / 3.0)) / pdur)
            tstims[j] = delay + ipi * j
        if stim['PT'] > 0.0:
            send = delay + ipi * j
            for i in range(send + posttest, len(w)):
                w[i] += stim['amp'] * (1.0 - exp(-i / (pdur / 3.0))) * exp(-(i - (pdur / 3.0)) / pdur)
    return (w, maxt, tstims)


def findPoint(x, y, peakindex, val, direction='left', limits=None):
    F = scipy.interpolate.UnivariateSpline(x, y, s=0) # declare function
    #To find x at y then do:
    istart = 0
    iend = len(y)
    if limits is not None:
        istart = int(limits[0] / limits[2])
        iend = int(limits[1] / limits[2])
    yToFind = y[peakindex] * val
    if direction == 'left':
        yreduced = numpy.array(y[istart:peakindex]) - yToFind
        try:
            Fr = scipy.interpolate.UnivariateSpline(x[istart:peakindex], yreduced, s=0)
        except:
            print 'findPoint: insufficient time points for analysis'
            print 'istart, peakindex: ', istart, peakindex
            print 'ytofine: ', yToFind
            res = float('nan')
            return (res)
        res = Fr.roots()
        if len(res) > 1:
            res = res[-1]
    else:
        yreduced = numpy.array(y[peakindex:iend]) - yToFind
        try:
            Fr = scipy.interpolate.UnivariateSpline(x[peakindex:iend], yreduced, s=0)
        except:
            print 'findPoint: insufficient time points for analysis'
            res = float('nan')
            return (res)
        res = Fr.roots()
        if len(res) > 1:
            res = res[0]
            #pdb.set_trace()
    try:
        res.pop()
    except:
        pass
    if not res: # tricky - an empty list is False, but does not evaluate to False
        res = float('nan') # replace with a NaN
    else:
        res = float(res) # make sure is just a simple number (no arrays)
    return (res)


def testSynapses(TargetCellName='bushy'):
    """ Test the synapse function. 
    v_pre is the presynaptic voltage
        v_soma is the postsynaptic voltage
        v_calyx is the calyx voltage
        coh is the calyx of held state
        isyn is the synaptic current
    """
    from scipy import interpolate
    import nrnlibrary.Cells as Cells

    h.nrn_load_dll(os.getcwd() + '../../PyNeuronLibrary/mechanisms/i386/libnrnmech.la')
    h.load_file("stdrun.hoc")
    vec = {}
    for var in ['v_pre', 'v_soma', 'i_soma', 'v_calyx', 'coh', 't', 'C0', 'C1', 'i_stim',
                'C2', 'C3', 'D', 'O1', 'O2', 'D1', 'D2', 'D3', 'Open', 'nmOpen', 'amOpen']:
        vec[var] = h.Vector()
    h.celsius = 34.0
    VCLAMP = True
    glyPlot = False
    releasePlot = True
    AN_Po_Ratio = 23.2917 # ration of open probabilities for AMPA and NMDAR's at peak currents
    AMPA_Max_Po = 0.44727
    NMDARatio = 0.0
    #
    # set parameters according to the target cell type
    #

    if TargetCellName == 'tstellate':
        TargetCell = Cells.tstellate_rothman(debug=True, ttx=False, nav11=True,
                                             species='mouse') # make a postsynaptic cell
        nANTerminals = 6
        nANTerminals_ReleaseZones = 1
        AN_gMax = 4600.0 / AMPA_Max_Po # correct because Po in model is not 1.0

    if TargetCellName == 'dstellate': # similar to t-stellate in abasence of other data
        TargetCell = Cells.dstellate_eager(debug=True, ttx=False) # make a postsynaptic cell
        nANTerminals = 12
        nANTerminals_ReleaseZones = 1
        AN_gMax = 4600.0 / AMPA_Max_Po # correct because Po in model is not 1.0

    if TargetCellName == 'bushy':
        nANTerminals = 3
        nANTerminals_ReleaseZones = 100
        # normally an_gmax would be about 22 nS (22000) total
        # However, if we do this as a "per release site" conductance, we have ~100 pA / 60 mV = 1.7 nS
        #AN_gMax = 22000.0/(AMPA_Max_Po) # correct because Po in model is not 1.0
        AN_gMax = 1700.0 / (AMPA_Max_Po)
        (TargetCell, [initseg, axn, internnode]) = Cells.bushy(debug=True, ttx=True) # make a postsynaptic cell

    #
    # make a multisite synapse (consists of an axon section and a coh point process)
    # synapse is like a calyx or end bulb of Held
    #
    calyx_zones = 1
    debug = False
    plotFocus = 'EPSC'
    axon = template_SGAxon()
    psdType = 'ampa'
    synType = 'epsc'
    ipscdynamics = 1

    ANTerminals_Delay = 0.0 # latency between AP and mode of release distribution, in milliseconds.    
    ANTerminals_Latency = 0.5 # latency 
    thresh = -30 # mV - AP detection on the presynaptic side.

    # set parameter control for the stochastic release of vesicles...
    # this structure is passed to stochastic synapses, and replaces several variables 
    # that were previously defined in the call to that function.
    vPars = Params(LN_Flag=1, LN_t0=10.0, LN_A0=0.05, LN_tau=35, LN_std=0.05,
                   Lat_Flag=1, Lat_t0=10.0, Lat_A0=0.140, Lat_tau=21.5,
                   delay=ANTerminals_Delay, latency=ANTerminals_Latency)

    printParams(vPars)
    if hasattr(TargetCell, 'len') and len(TargetCell) > 1: # handle possible multi-segment cells
        TC = TargetCell[0] # soma
        SOMA = TargetCell[0]
    else:
        TC = TargetCell
        SOMA = TargetCell
    if TargetCellName == 'dstellate': # well, these always have a dendrite... 
        TC = TargetCell[0]
        SOMA = TargetCell

    for i in range(0, nANTerminals_ReleaseZones): # make a synapse monitor for each release zone
        vec['isyn%03d' % i] = h.Vector(nANTerminals_ReleaseZones, 1000)
        #
    # stochastic_synapses builds the multisite synapse.     
    #
    (calyx, coh, psd, cleft, nc2, par) = stochastic_synapses(h, axon=axon,
                                                             targetcell=TC, cellname=TargetCellName,
                                                             nFibers=nANTerminals, nRZones=nANTerminals_ReleaseZones,
                                                             eRev=0, debug=False,
                                                             thresh=thresh, psdtype=psdType, gmax=AN_gMax, gvar=0.3,
                                                             NMDARatio=0.0, Identifier=1,
                                                             stochasticPars=vPars) # set gVar to 0 for testing
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

    # adjust NMDA receptors
    #print dir(axon)
    k = 0
    kNMDA = -1
    kAMPA = -1
    #    for p in psd:
    #        print dir(p)
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
    (secmd, maxt, tstims) = make_pulse(stim, pulsetype='square')
    istim.delay = 0
    istim.dur = 1e9 # these actually do not matter...
    istim.iMax = 0.0
    # istim current pulse train
    vec['i_stim'] = h.Vector(secmd)

    # set up recordings

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

    h.tstop = 200.0 # duration of a run
    tvec = arange(0, h.tstop, h.dt)
    h.init() # set up the parameters
    for nrep in xrange(1): # could do multiple runs.... 
        h.run()
        k = 0
        if nrep is 0: # save the soma current
            isoma = numpy.zeros_like(vec['isyn000'].to_python())
        for p in psd:
            thissyn = 'isyn%03d' % k
            if thissyn in vec.keys():
                isoma = isoma + vec[thissyn].to_python()
                k = k + 1
    print 'Synapse.py: all runs done'
    nreq = 0
    nrel = 0
    ntrel = numpy.zeros(nANTerminals)
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
        nmOmax = numpy.asarray(vec['nmOpen']).max()
        amOmax = numpy.asarray(vec['amOpen']).max()
        print 'Synapse.py: Max NMDAR Open Prob: %f   AMPA Open Prob: %f\n' % (nmOmax, amOmax)
        nmImax = numpy.asarray(vec['isyn%03d' % kNMDA]).max()
        amImax = numpy.asarray(vec['isyn%03d' % kAMPA]).max()
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
    tpsc = numpy.arange(0, ipi + textend, h.dt)
    ipsc = numpy.zeros((stim['NP'], pscpts))
    mpl.figure(num=220, facecolor='w')
    gpsc = mpl.subplot2grid((5, 2), (0, 0), rowspan=2, colspan=2)
    ti = numpy.asarray(vec['t'])
    psc_20_lat = numpy.zeros((stim['NP'], 1)) # latency to 20% of rising amplitude
    psc_80_lat = numpy.zeros((stim['NP'], 1)) # latency to 80% of rising amplitude
    psc_hw = numpy.zeros((stim['NP'], 1)) # width at half-height
    psc_rt = numpy.zeros((stim['NP'], 1)) # 20-80 rise time
    tp = numpy.zeros((stim['NP'], 1))
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
        psc_pk = numpy.argmax(ipsc[i, minStart:]) # position of the peak
        psc_pk = psc_pk + minStart - 1
        print 'i, pscpk, ipsc[i,pscpk]: ', i, psc_pk, ipsc[i, psc_pk]
        #       print 'minLat: %f   ipi+textend: %f, hdt: %f' % ((minLat, ipi+textend, h.dt))
        if psc_pk == 0:
            continue
        pkval = ipsc[i, psc_pk]
        psc_20_lat[i] = findPoint(tpsc, ipsc[i, :], psc_pk, 0.2, direction='left', limits=(minLat, ipi + textend, h.dt))
        psc_80_lat[i] = findPoint(tpsc, ipsc[i, :], psc_pk, 0.8, direction='left', limits=(minLat, ipi + textend, h.dt))
        psc_50l = findPoint(tpsc, ipsc[i, :], psc_pk, 0.5, direction='left', limits=(minLat, ipi + textend, h.dt))
        psc_50r = findPoint(tpsc, ipsc[i, :], psc_pk, 0.5, direction='right', limits=(minLat, ipi + textend, h.dt))
        if not math.isnan(psc_20_lat[i]) and not math.isnan(psc_80_lat[i]):
            psc_rt[i] = psc_80_lat[i] - psc_20_lat[i]
        else:
            psc_rt[i] = numpy.nan
        if not math.isnan(psc_50r) and not math.isnan(psc_50l):
            psc_hw[i] = float(psc_50r) - float(psc_50l)
            gpsc.plot(psc_50l, pkval * 0.5, 'k+')
            gpsc.plot(psc_50r, pkval * 0.5, 'k+')
            gpsc.plot(tpsc, ipsc[i, :].T)
        else:
            psc_hw[i] = numpy.nan
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
    RT_mean2080_early = numpy.ma.masked_invalid(psc_rt[analysisWindow[0]]).mean()
    RT_mean2080_late = numpy.ma.masked_invalid(psc_rt[analysisWindow[1]]).mean()
    Lat_mean20_early = numpy.ma.masked_invalid(psc_20_lat[analysisWindow[0]]).mean()
    Lat_mean20_late = numpy.ma.masked_invalid(psc_20_lat[analysisWindow[1]]).mean()
    HW_mean_early = numpy.ma.masked_invalid(psc_hw[analysisWindow[0]]).mean()
    HW_mean_late = numpy.ma.masked_invalid(psc_hw[analysisWindow[1]]).mean()
    print "Means: --------------"
    print RT_mean2080_early
    print Lat_mean20_early
    print HW_mean_early
    print 'Early:   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_mean2080_early, Lat_mean20_early,
                                                                              HW_mean_early)
    print 'Late :   RT {0:7.3f} ms   Lat {1:7.3f} ms   HW {2:7.3f} ms'.format(RT_mean2080_late, Lat_mean20_late,
                                                                              HW_mean_late)
    RT_std2080_early = numpy.ma.masked_invalid(psc_rt[analysisWindow[0]]).std()
    RT_std2080_late = numpy.ma.masked_invalid(psc_rt[analysisWindow[1]]).std()
    Lat_std20_early = numpy.ma.masked_invalid(psc_20_lat[analysisWindow[0]]).std()
    Lat_std20_late = numpy.ma.masked_invalid(psc_20_lat[analysisWindow[1]]).std()
    HW_std_early = numpy.ma.masked_invalid(psc_hw[analysisWindow[0]]).std()
    HW_std_late = numpy.ma.masked_invalid(psc_hw[analysisWindow[1]]).std()
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
            dis = numpy.array(coh[j].EventDist)[0:nev] # actual events, not original distribution
            tim = numpy.array(coh[j].EventTime)[0:nev]
            #            print "dis[%03d]: nev=%d " % (j, nev),
            #            print "# of spike Requests detected: %d \n" % (coh[j].nRequests)
            #print dis
            (hist, binedges) = numpy.histogram(dis, 25, range=(0.0, 5.0))
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

################################################################################
#
# If we call this directly, provide a test of the synapse functions
# - see below to switch cells
#
################################################################################

if __name__ == "__main__":
    testSynapses()