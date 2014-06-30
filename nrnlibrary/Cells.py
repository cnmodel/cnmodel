#!/usr/bin/python
#
# Cell definitions for models.
#
# This file includes a number of different cell definitions and default
# conductances for point models. Most are models from the lab for neurons
#of the cochlear nucleus.
# Pyramidal cell from the DCN: pyr (Kanold and Manis, 1999, 2001, 2005)
# Bushy cell from the VCN: bushy (Rothman and Manis, 2003abc) Type II
# T-stellate cell from teh VCN: tstellate (Rothman and Manis, 2003abc) Type I
#
# Paul B. Manis, Ph.D. 2009 (August - November 2009)
#

import os
import os.path
from neuron import h
from neuron import *
from math import sqrt
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

runQuiet = True


def getcelltypes():
    """
    Return the types of cells that area generally avaialble.
    This is sorely outdated for the files here.
    """
    return({'pyramidal': (pyr, 'pyrprc.txt',
                'PYR-PRC.py output file for pyramidal'),
            'cartwheel': (cartwheel, 'cwprc.txt',
                'CW-PRC.py output file for Cartwheel'),
            'bushy': (bushy, 'bushyprc.txt',
                'BU-PRC.py output file for Bushy'),
            'tstellate': (tstellate_rothman, 'tstellateprc.txt',
                'ST-PRC.py output file for T-Stellate')})
# just return dictionary of the cells defined in this file

#-----------------------------------------------------------------------------
# Standard Hodkgin-Huxley model, small cell


def hh(debug=False, message=None):
    """
    Standard Hodgkin-Huxley mechanisms from NEURON
    """
    soma = h.Section() # one compartment of about 29000 um2
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential
    c_m = 1.0
    scalefactor = 1.0 # This determines the relative size of the cell
    rinsf = 1.0           # input resistance adjustment (also current...)
    totcap = 20.0 # scalefactor * 1.0 # cap in pF for cell
    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / c_m # pf -> uF, cm = 1uf/cm^2 nominal
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd

    seg = soma
    seg.insert('hh')
    seg.insert('pas')
    if not runQuiet:
        if message is None:
            print "<< Standard HH model created >>"
        else:
            print message
    return(soma, (None, None, None))


def hasenstaub(debug=False, ttx=False, message=None, pump=False):
    v_potassium = -70       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential
    v_chloride = -20
    cm = 1.0
    Ra = 150
    scalefactor = 1.0 # This determines the relative size of the cell
    rinsf = 1.0           # input resistance adjustment (also current...)
    totcap = 1e5 # in pF

    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    soma = h.Section() # one compartment of about 29000 um2
    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd
    print 'soma: (um len, um dia, area in cm2)', soma.L, soma.diam, somaarea

    soma.insert('hstb')
    gkb = nstomho(7000, somaarea) * scalefactor
    gnab = nstomho(1500, somaarea) * scalefactor
    gl = nstomho(100, somaarea) * scalefactor
    soma().hstb.gbar = gnab
    soma().hstb.gbar = gkb
    soma().hstb.gl = gl
    print 'gkb, gnab: ', gkb, gnab
    if pump:
        soma.insert('k_conc')

        ki0_k_ion = 140
        soma().ki = ki0_k_ion
        soma().ki0_k_conc = ki0_k_ion
        soma().beta_k_conc = 0.075

        soma.insert('na_conc')
        nai0_na_ion = 5
        soma().nai = nai0_na_ion
        soma().nai0_na_conc = nai0_na_ion
        soma().beta_na_conc = 0.075

        soma.insert('nakpump')
        soma().nakpump.inakmax = 8
        soma().nao = 145
        soma().ko = 5
        soma().nakpump.Nai_inf = 5
        soma().nakpump.Ki_inf = 140
        soma().nakpump.ATPi = 5
    soma.ek = v_potassium
    soma().v = -60.
    if ttx:
        gbar = 0.0
        soma().na.gbar = gbar
        soma.ena = 50

    if not runQuiet:
        if message is None:
            print "<< bushy: JSR bushy cell model created >>"
        else:
            print message
    return(soma, [None, None, None])




def pyr(scalefactor=1.0, debug=False):

    """"""
    ndend = 1
    soma = h.Section()
    dend = h.Section()

    cm = 1.0
    scalefactor = 1.0 # This determines the relative size of the cell
    rinsf = 1.0           # input resistance adjustment (also current...)
    totcap = scalefactor * 12.0 # cap in pF for cell
    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    gnab = nstomho(350, somaarea) * scalefactor
    gnap = 0.0
    gkb = nstomho(80, somaarea) * scalefactor # used to be 20?
    gkfb = nstomho(150, somaarea) * scalefactor
    gksb = nstomho(40, somaarea) * scalefactor
    glb = nstomho(2.8, somaarea) * scalefactor
    ghb = nstomho(3, somaarea) * scalefactor
    gkpksk = nstomho(0, somaarea) * scalefactor
    gkir = nstomho(0, somaarea) * scalefactor # incude KIR here, but set to 0
    stdrin = 300 * rinsf / scalefactor
    stdrmp = -60
    if debug:
        print ("in Mhos: gna: %9.3g  gkb: %9.3g  gkab: %9.3g  gksb: %9.3g"
            % (gnab, gkb, gkfb, gksb))
        print ("glb: %9.3g  ghb: %9.3g gkpksk: %9.3g" %
            (glb, ghb, gkpksk))

# set up soma like a pyramidal cell
    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd # these are expressed in microns...
    soma.insert('pyr')
    soma.insert('kpksk')
    soma.insert('cadiff') # diffusion
    soma.insert('cap') # p-type calcium current
    #soma.insert('nacum') # sodium accumulation (yes!)
    soma.insert('nakpump') # and of course a pump to handle it.
    soma.insert('k_conc')
    soma.insert('na_conc')
    # soma.insert('kna')

    seg = soma()
    seg.kpksk.gbar = gkpksk
    seg.cap.pcabar = 0.00002
    seg.pyr.gbar = gnab
    seg.pyr.gbar = gnap
    seg.pyr.gbar = gkb
    seg.pyr.gbar = gkfb
    seg.pyr.gbar = gksb
    seg.pyr.gl = glb
    seg.pyr.gbar = ghb
   # seg.pyr.gbar = gkir
    seg.pyr.kif_ivh = -89.6
    if not runQuiet:
        print " "
        print "<< PYR: POK Pyramidal Cell created >>"
        print " "
    return(soma) # return the section we defined


def cartwheel(debug=False):
    soma = h.Section() # one compartment of about 29000 um2
    soma.nseg = 1
    soma.diam = 96
    soma.L = 96
    cm = 1
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential

    seg = soma
    seg.insert('naRsg')
    seg.insert('kpkj')
    seg.insert('kpkj2')
    seg.insert('kpkjslow')
    seg.insert('bkpkj')
    seg.insert('kpksk')
    seg.insert('cadiff')
    seg.insert('cap')
    seg.insert('lkpkj')
    seg.insert('hpkj')
    seg.ena = 60
    seg.ek = -80
    s = soma()
    s.kpksk.gbar = 0.002
    if not runQuiet:
        print "<< cartwheel: Raman Purkinje cell model (modified) created >>"
    return(soma)


def sgc(debugFlag=False, ttx=False, message=None, nach='jsrna',
        species='mouse', chlist=None):
    """
    Definitions for spiral ganglion cells
    """
    if chlist is None:
        chlist = ['ih', 'klt', 'kht', 'na']
    v_potassium = -84       # potassium reversal potential
    v_sodium = 55           # sodium reversal potential
    v_chloride = -20
    v_eh = -41.3 # from R. Davis papers
    cm = 1.0
    scalefactor = 1.0 # This determines the relative size of the cell
    rinsf = 1.0           # input resistance adjustment (also current...)
    if species == 'guineapig':
        totcap = scalefactor * 12.0 # cap in pF for cell
    if species == 'mouse':
        totcap = scalefactor * 12.0 # cap in pF for cell from JS-S
    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    soma = h.Section() # one compartment of about 29000 um2
    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd

    if 'klt' in chlist:
        soma.insert('klt')
    if 'kht' in chlist:
        soma.insert('kht')
    if 'na' in chlist:
        if nach == 'jsrna':
            soma.insert('jsrna')
        elif nach == 'nav11':
            soma.insert('nav11')
        else:
            soma.insert('na')
    #soma.insert('ka')
    if 'ih' in chlist:
        soma.insert('ihsgc')
    #soma.insert('hcno')
    soma.insert('leak')
    if 'klt' in chlist or 'klt' in chlist:
        soma.ek = v_potassium

    s = soma()
    gnabar = nstomho(1000.0, somaarea) * scalefactor
    if 'na' in chlist:
        if ttx:
            gnabar = 0.0
        if nach == 'jsrna':
            s.jsrna.gna = gnabar
            soma.ena = 50
        elif nach == 'nav11':
            s.nav11.gnatbar = gnabar * 0.5
            soma.ena = 50
            s.nav11.vsna = 4.3
            print "sgc using inva11"
        else:
            s.na.gbar = gnabar
            soma.ena = 50
    if species == 'mouse':
        # use conductance levels from Cao et al.,  J. Neurophys., 2007.
        #if debug:
        print 'Mouse sgc cell'
        #s.nav11.gbar = gnabar*0.3
        if 'kht' in chlist:
            s.kht.gbar = nstomho(58.0, somaarea) * scalefactor
        if 'klt' in chlist:
            s.klt.gbar = nstomho(80.0, somaarea) * scalefactor
            # nstomho(200.0, somaarea) * scalefactor
        #s.ka.gkabar = nstomho(0.0, somaarea) * scalefactor
        if 'ih' in chlist:
            s.ihsgc.gbar = nstomho(10.0, somaarea) * scalefactor
            s.ihsgc.eh = v_eh
        #s.hcno.gbar = 0.0
        s.lek.gbar = nstomho(0.5, somaarea) * scalefactor
        vm0 = -63.6
    if species == 'guineapig-sgc-II':
         # guinea pig data from Rothman and Manis, 2003, type II
        if debugFlag:
            print 'Guinea pig sgc cell type II'
        if 'kht' in chlist:
            s.kht.gbar = nstomho(150.0, somaarea) * scalefactor
        if 'klt' in chlist:
            s.klt.gbar = nstomho(200.0, somaarea) * scalefactor
            # nstomho(200.0, somaarea) * scalefactor
        #s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
        if 'ih' in chlist:
            s.ihsgc.gbar = nstomho(1.0, somaarea) * scalefactor
        #s.hcno.gbar = 0.0
        s.lek.gbar = nstomho(2.0, somaarea) * scalefactor
        vm0 = -63.6

    if not runQuiet:
        if message is None:
            print "<< sgc: SGC cell model created >>"
        else:
            print message
    return(soma, [None, None, None])


def bushy(debug=False, ttx=False, message=None, nach='jsrna',
    species='mouse', axon=False, dendrite=False,
    newModFiles=False, pump=False):
    v_potassium = -70       # potassium reversal potential
#    v_sodium = 55           # sodium reversal potential
#    v_chloride = -20
    c_m = 1.0
    R_a = 150
    scalefactor = 1.0 # This determines the relative size of the cell
#    rinsf = 1.0           # input resistance adjustment (also current...)
    if species == 'guineapig-bushy-II' or species == 'guineapig-bushy-II-I':
        totcap = scalefactor * 12.0 # cap in pF for cell
#        axonsf = 0.57
#        isdist = 0
    if species == 'mouse':
        totcap = scalefactor * 26.0 # cap in pF for cell from Cao et al., 2007
#        axonsf = 0.57
#        isdist = 0
    if species == 'cat':
        totcap = scalefactor * 35.0
#        axonsf = 1.0
#        isdist = 0

#    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / c_m # pf -> uF, cm = 1uf/cm^2 nominal
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    soma = h.Section() # one compartment of about 29000 um2
    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd

    if newModFiles:
        soma.insert('klt2')
        soma.insert('kht2')

    if pump:
        soma.insert('k_conc')

        ki0_k_ion = 140
        soma().ki = ki0_k_ion
        soma().ki0_k_conc = ki0_k_ion
        soma().beta_k_conc = 0.075

        soma.insert('na_conc')
        nai0_na_ion = 5
        soma().nai = nai0_na_ion
        soma().nai0_na_conc = nai0_na_ion
        soma().beta_na_conc = 0.075

        soma.insert('nakpump')
        soma().nakpump.inakmax = 8
        soma().nao = 145
        soma().ko = 5
        soma().nakpump.Nai_inf = 5
        soma().nakpump.Ki_inf = 140
        soma().nakpump.ATPi = 5
    else:
        soma.insert('klt')
        soma.insert('kht')

    if nach == 'jsrna':
        soma.insert('jsrna')
    elif nach == 'nav11':
        soma.insert('nav11')
    else:
        soma.insert('na')
    #soma.insert('ka')
    soma.insert('ihvcn')
    #soma.insert('hcno')
    soma.insert('leak')
    soma.ek = v_potassium
    gnabar = nstomho(1000.0, somaarea) * scalefactor
    if ttx:
        gnabar = 0.0
    if nach == 'jsrna':
        soma().jsrna.gna = gnabar
        soma.ena = 50
    elif nach == 'nav11':
        soma().nav11.gnatbar = gnabar * 0.5
        soma.ena = 50
        soma().nav11.vsna = 4.3
        print "bushy using inva11"
    else:
        soma().na.gnabar = gnabar
        soma.ena = 50

    vm0 = bushy_species_scaling(soma, species, newModFiles, somaarea, scalefactor, debug)

    maindend = None
    secdend = None
    if dendrite:
        print 'Adding dendrite to Bushy model'
        section = h.Section
        maindend = section(cell=soma)
        maindend.connect(soma)
        maindend.nseg = 10
        maindend.L = 100.0
        maindend.diam = 2.5
        maindend.insert('klt')
        maindend.insert('ihvcn')
        if newModFiles:
            maindend().klt2.gbar = soma().klt2.gbar / 2.0
        else:
            maindend().klt.gbar = soma().klt.gbar / 2.0

        maindend().ihvcn.gbar = soma().ihvcn.gbar / 2.0

        maindend.cm = c_m
        maindend.Ra = R_a
        nsecd = range(0, 5)
        secdend = []
        for ibd in nsecd:
            secdend.append(section(cell=soma))
        for ibd in nsecd:
            secdend[ibd].connect(maindend)
            secdend[ibd].diam = 1.0
            secdend[ibd].L = 15.0
            secdend[ibd].cm = c_m
            secdend[ibd].Ra = R_a
        #h.topology()

    if not runQuiet:
        if message is None:
            print "<< bushy: JSR bushy cell model created >>"
        else:
            print message
    return(soma, [maindend, secdend, None])

def bushy_species_scaling(soma, species, newModFiles, somaarea, scalefactor, debug):
    if species == 'mouse':
        # use conductance levels from Cao et al.,  J. Neurophys., 2007.
        #if debug:
        print 'Mouse bushy cell'
        #s.nav11.gbar = gnabar*0.3
        if newModFiles:
            soma().kht2.gbar = nstomho(58.0, somaarea) * scalefactor
            soma().klt2.gbar = nstomho(80.0, somaarea) * scalefactor
        else:
            soma().kht.gbar = nstomho(58.0, somaarea) * scalefactor
            soma().klt.gbar = nstomho(80.0, somaarea) * scalefactor
        #s.ka.gkabar = nstomho(0.0, somaarea) * scalefactor
        soma().ihvcn.gbar = nstomho(30.0, somaarea) * scalefactor
        #s.hcno.gbar = 0.0
        soma().leak.gbar = nstomho(2.0, somaarea) * scalefactor
        vm0 = -63.6
    if species == 'guineapig-bushy-II':
        # guinea pig data from Rothman and Manis, 2003, type II
        if debug:
            print 'Guinea pig bushy cell type II'
        if newModFiles:
            soma().kht2.gbar = nstomho(150.0, somaarea) * scalefactor
            soma().klt2.gbar = nstomho(200.0, somaarea) * scalefactor
        else:
            soma().kht.gbar = nstomho(150.0, somaarea) * scalefactor
            soma().klt.gbar = nstomho(200.0, somaarea) * scalefactor
        #s.ka.gkabar = nstomho(0.0, somaarea) * scalefactor
        soma().ihvcn.gbar = nstomho(20.0, somaarea) * scalefactor
        #s.hcno.gbar = 0.0
        soma().lek.gbar = nstomho(2.0, somaarea) * scalefactor
        vm0 = -63.6
    if species == 'guineapig-bushy-II-I':
         # guinea pig data from Rothman and Manis, 2003, type II=I
        if debug:
            print 'Guinea pig bushy cell type II-I'
        if newModFiles:
            soma().kht2.gbar = nstomho(150.0, somaarea) * scalefactor
            soma().klt2.gbar = nstomho(35.0, somaarea) * scalefactor
        else:
            soma().kht.gbar = nstomho(150.0, somaarea) * scalefactor
            soma().klt.gbar = nstomho(35.0, somaarea) * scalefactor
        #s.ka.gkabar = nstomho(0.0, somaarea) * scalefactor
        soma().ihvcn.gbar = nstomho(3.5, somaarea) * scalefactor
        #s.hcno.gbar = 0.0
        soma().lek.gbar = nstomho(2.0, somaarea) * scalefactor
        vm0 = -63.8
    return vm0


# def bushy_addAxon(debug=False, ttx=False, message=None, nach=None):
# 
#     nnodes = range(5)
# #    print nnodes
#     axnode = []
#     internode = []
#     section = h.Section
#     initsegment = section(cell=soma)
#     initsegment.connect(soma)
#     #for i in range(nnodes):
#     #    axnode.append(Section(cell=soma))
#     #    internode.append(Section(cell=soma))
#     for i in nnodes:
#         axnode.append(Section(cell=soma))
#         internode.append(Section(cell=soma))
# 
#     axnode[0].connect(initsegment)
#     for i in nnodes:
#         internode[i].connect(axnode[i])
#         if i < nnodes[-1]:
#             axnode[i + 1].connect(internode[i])
# # create an initial segment
#     ninitseg = 21
#     initsegment.nseg = ninitseg
#     initsegment.diam = 4.0 * axonsf
#     initsegment.L = 36.0 * axonsf
#     initsegment.cm = cm
#     initsegment.Ra = Ra
#     initsegment.insert('na')
#     initsegment.insert('kht')
#     initsegment.insert('klt')
#     initsegment.insert('ihvcn')
#     initsegment.insert('leak')
#     gnamax = nstomho(6000.0, somaarea) * scalefactor
#     gnamin = 0.0 * gnamax
# 
#     gnastep = (gnamax - gnamin) / ninitseg
#     ip = 0
#     for inseg in initsegment:
#         ina = gnamin + ip * gnastep
#         print 'seg %d ina = %9.6f' % (ip, ina)
#         inseg.na.gbar = ina
#         inseg.klt.gbar = 0.2 * nstomho(200.0, somaarea) * scalefactor
#         inseg.kht.gbar = nstomho(150.0, somaarea) * scalefactor
#         inseg.ihvcn.gbar = 0.0 * nstomho(20.0, somaarea) * scalefactor
#         inseg.lek.gbar = nstomho(2.0, somaarea) * scalefactor
#         inseg.ena = v_sodium
#         inseg.ek = v_potassium
#         ip = ip + 1
# 
#     for i in nnodes:
#         axnode[i] = loadaxnodes(axnode[i], somaarea, scalefactor)
#         internode[i] = loadinternodes(internode[i], somaarea, scalefactor)
# 
#    # h.topology()
#     vm0 = -63.6
#     if not runQuiet:
#         if message is None:
#             print "<< bushy: JSR bushy cell model with axon created >>"
#         else:
#             print message
#     return(soma, [initsegment, axnode, internode])


def bushy_waxon(debug=False, ttx=False, message=None, nach=None):
    v_potassium = -70       # potassium reversal potential
    v_sodium = 55           # sodium reversal potential
    Ra = 150
    cm = 1.0

    species = 'mouse'
    if species is 'cat':
        axonsf = 1.0
        isdist = 0
    if species is 'mouse':
        axonsf = 0.57
        isdist = 0
    scalefactor = axonsf # This determines the relative size of the cell
#    rinsf = 1.0           # input resistance adjustment (also current...)
    totcap = scalefactor * 12.0 # cap in pF for cell
    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    nnodes = range(5)
#    print nnodes
    axnode = []
    internode = []
    Section = h.Section
    soma = Section() # one compartment of about 29000 um2
    initsegment = Section(cell=soma)
    initsegment.connect(soma)
    #for i in range(nnodes):
    #    axnode.append(Section(cell=soma))
    #    internode.append(Section(cell=soma))
    for i in nnodes:
        axnode.append(Section(cell=soma))
        internode.append(Section(cell=soma))

    axnode[0].connect(initsegment)
    for i in nnodes:
        internode[i].connect(axnode[i])
        if i < nnodes[-1]:
            axnode[i + 1].connect(internode[i])

    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd
    soma.Ra = Ra
    soma.cm = cm
    soma.insert('klt')
    soma.insert('kht')
    soma.insert('na')
    soma.insert('ihvcn')
    soma.insert('leak')
    soma.ena = v_sodium
    soma.ek = v_potassium
#    print dir(soma)
    for s in soma:
        if ttx is False:
            s.na.gbar = 0. * nstomho(1000.0, somaarea) * scalefactor
        else:
            s.na.gbar = 0.0
        s.kht.gbar = nstomho(150.0, somaarea) * scalefactor
        s.klt.gbar = nstomho(200.0, somaarea) * scalefactor
        s.ihvcn.gbar = nstomho(20.0, somaarea) * scalefactor
        s.lek.gbar = nstomho(2.0, somaarea) * scalefactor

# create an initial segment
    ninitseg = 21
    initsegment.nseg = ninitseg
    initsegment.diam = 4.0 * axonsf
    initsegment.L = 36.0 * axonsf
    initsegment.cm = cm
    initsegment.Ra = Ra
    initsegment.insert('na')
    initsegment.insert('kht')
    initsegment.insert('klt')
    initsegment.insert('ihvcn')
    initsegment.insert('leak')
    gnamax = nstomho(6000.0, somaarea) * scalefactor
    gnamin = 0.0 * gnamax

    gnastep = (gnamax - gnamin) / ninitseg
    ip = 0
    for inseg in initsegment:
        ina = gnamin + ip * gnastep
        print 'seg %d ina = %9.6f' % (ip, ina)
        inseg.na.gbar = ina
        inseg.klt.gbar = 0.2 * nstomho(200.0, somaarea) * scalefactor
        inseg.kht.gbar = nstomho(150.0, somaarea) * scalefactor
        inseg.ihvcn.gbar = 0.0 * nstomho(20.0, somaarea) * scalefactor
        inseg.lek.gbar = nstomho(2.0, somaarea) * scalefactor
        inseg.ena = v_sodium
        inseg.ek = v_potassium
        ip = ip + 1

    for i in nnodes:
        axnode[i] = loadaxnodes(axnode[i], somaarea, scalefactor)
        internode[i] = loadinternodes(internode[i], somaarea, scalefactor)

   # h.topology()
    vm0 = -63.6
    if not runQuiet:
        if message is None:
            print "<< bushy: JSR bushy cell model with axon created >>"
        else:
            print message
    return(soma, [initsegment, axnode, internode])


def loadaxnodes(axnode, somaarea, scalefactor, nodeLength=2.5,
    nodeDiameter=2.0):
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential
    Ra = 150
    cm = 1.0
    axnode.nseg = 1
    axnode.L = nodeLength * scalefactor
    axnode.diam = nodeDiameter * scalefactor
    axnode.Ra = Ra
    axnode.cm = cm
    axnode.insert('na')
    axnode.insert('kht')
    axnode.insert('klt')
    axnode.insert('leak')
    axnode.insert('ihvcn')
    for ax in axnode:
    #    print "   segment: ", ax
        ax.na.gbar = nstomho(1000.0, somaarea) * scalefactor
        ax.kht.gbar = nstomho(150.0, somaarea) * scalefactor
        ax.klt.gbar = nstomho(200.0, somaarea) * scalefactor
        ax.ihvcn.gbar = 0
        ax.lek.gbar = nstomho(2.0, somaarea) * scalefactor
        ax.ena = v_sodium
        ax.ek = v_potassium
    return(axnode)


def loadinternodes(internode, somaarea, scalefactor, internodeLength=1000,
    internodeDiameter=10):
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential
    Ra = 150
    cm = 1.0

    internode.nseg = 20
    internode.L = internodeLength * scalefactor
    internode.diam = internodeDiameter * scalefactor
    internode.Ra = Ra
    internode.cm = 0.002
    internode.insert('na')
    internode.insert('kht')
    internode.insert('leak')
    for inno in internode:
        inno.lek.gbar = nstomho(0.002, somaarea) * scalefactor
        inno.na.gbar = 0 * nstomho(500.0, somaarea) * scalefactor
        inno.kht.gbar = 0 * nstomho(150.0, somaarea) * scalefactor
        inno.ek = v_potassium
        inno.ena = v_sodium
        inno.leak.e = -80
    return(internode)


def tstellate_rothman(debug=False, ttx=False, message=None,
    species='guinea pig', nav11=False):
    soma = h.Section() # one compartment of about 29000 um2
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential

    cm = 1.0
    scalefactor = 1.0 # This determines the relative size of the cell
    rinsf = 1.0           # input resistance adjustment (also current...)
    if species == 'guinea pig':
        totcap = scalefactor * 12.0 # cap in pF for cell
    elif species == 'mouse':
        totcap = scalefactor * 25.0 # adjusted...
    else:
        raise ValueError('Unrecognized species in call to tstellate_rothman')

    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd

    seg = soma
    seg.insert("kht")
    if nav11 is False:
        seg.insert('na')
        seg.ena = 50
    else:
        seg.insert('nav11')
        seg.ena = 50
    seg.insert('ka')
    seg.insert('ihvcn')
    seg.insert('leak')
    seg.ek = -84
    s = soma()
    if species == 'guinea pig':
        gnamax = 1000.0
        s.kht.gbar = nstomho(150.0, somaarea) * scalefactor
        s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
        s.ihvcn.gbar = nstomho(0.5, somaarea) * scalefactor
        s.ihvcn.eh = -43 # Rodrigues and Oertel, 2006
        #print 'ih vcn vh: %f ' % (s.ihvcn.vh)
        s.lek.gbar = nstomho(2.0, somaarea) * scalefactor
        vm0 = -63.9
    elif species == 'mouse':
        gnamax = 800.0
        s.kht.gbar = nstomho(250.0, somaarea) * scalefactor
        s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
        s.ihvcn.gbar = nstomho(18.0, somaarea) * scalefactor
        s.ihvcn.eh = -43 # Rodrigues and Oertel, 2006
        s.lek.gbar = nstomho(8.0, somaarea) * scalefactor
        # yields input resistance of 74.2 Mohm measured from -60 to -70 mV
        vm0 = -60.0
    else:
        raise ValueError('species not recognized in tstellate_rothman')

    if not ttx:
        if nav11 is False:
            s.na.gna = nstomho(gnamax, somaarea) * scalefactor
        else:
            s.nav11.gbar = nstomho(1800.0, somaarea) * scalefactor
            s.nav11.vsna = 4.3 # was 8
    else:
        if nav11 is False:
            s.na.gna = 0.0
        else:
            s.nav11.gbar = 0.0

    if not runQuiet:
        if message is None:
            print "<< T-stellate: JSR Stellate Type 1 cell model created >>"
        else:
            print message
    return(soma)


def tstellate_rothman_nav11(debug=False, ttx=False, cs = False, message=None, dend=False):
    """
    T-stellate cell setup from Rothman and Manis, 2003, using nav11 sodium channel model
    ttx: if True, turns off sodium channels
    cs: if True, turns off K channels (e.g., cesium in pipette). 
    dend: if True, adds dendrites to the model, based roughly on White et al.,
    1994)
    NOTE: This has been modified from it's original from
    for use in simulating MOUSE stellate cells.
    """
    print ("T-STELLATE ROTHMAN",
        "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    soma = h.Section() # one compartment of about 29000 um2
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential

    cm = 1.0
    scalefactor = 1.0 # This determines the relative size of the cell
    rinsf = 1.0           # input resistance adjustment (also current...)
    totcap = scalefactor * 25.0 # cap in pF for cell
    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd

    dendrites=[]
    if dend is True:
        nDend = range(4) # these will be simple, unbranced, N=4 dendrites

        
#    print nnodes
        for i in nDend:
            dendrites.append(h.Section(cell=soma))
        for i in nDend:
            dendrites[i].connect(soma)
            dendrites[i].L = 200 # length of the dendrite (not tapered)
            dendrites[i].diam = 1.5 # dendritic diameter
            dendrites[i].nseg = 21 # # segments in dendrites
            dendrites[i].Ra = 150 # ohm.cm
            d = dendrites[i]
            ds = d()
            d.insert('kht')
            if cs is False:
                ds.kht.gbar = 0.005 # a little Ht
            else:
                ds.kht.gbar = 0.0
            d.insert('leak') # leak
            ds.lek.gbar = 0.0001
            d.insert('ihvcn') # some H current
            ds.ihvcn.gbar = 0.# 0.001
            ds.ihvcn.eh = -43.0
    seg = soma
    seg.insert('kht')
    seg.insert('nav11')
    seg.insert('ka')
    seg.insert('ihvcn')
    seg.insert('leak')
    seg.ena = 10
    seg.ek = -84
    s = soma()
    if ttx is False:
        s.nav11.gbar = nstomho(1800.0, somaarea) * scalefactor
    else:
        s.nav11.gbar = 0.0 # print s.nav11.gnat
    s.nav11.vsna = 4.3 # was 8

    if cs is False:
        s.kht.gbar = nstomho(200.0, somaarea) * scalefactor
        s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
    else:
        s.kht.gbar = 0.
        s.ka.gbar = 0.
    s.ihvcn.gbar = nstomho(18.0, somaarea) * scalefactor # was 10
    s.ihvcn.eh = -43 # Rodrigues and Oertel, 2006
#    print 'ih vcn vh: %f ' % (s.ihvcn.vh)
    s.lek.gbar = nstomho(2.0, somaarea) * scalefactor
    vm0 = -63.9
    if not runQuiet:
        if message is None:
            print ("<< T-stellate: JSR Stellate Type 1 cell model created",
            " - modified for mouse >>")
        else:
            print message
   # print dendrites
    return(soma, dendrites, None)


def tstellate_f(debug=False, ttx=False, message=None, dend=False):
    """ t-stellate model based on Rothman and Manis 2003, but with fast sodium
    channel """
    soma = h.Section() # one compartment of about 29000 um2
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential

    cm = 1.0
    scalefactor = 1.0 # This determines the relative size of the cell
    rinsf = 1.0           # input resistance adjustment (also current...)
    totcap = scalefactor * 20.0 # cap in pF for cell
    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd

    seg = soma
    seg.insert('kht')
    seg.insert('nav11')
    seg.insert('ka')
    # 'it' is not part of canonical model; 
    # just trying it to reproduce some data.
    #seg.insert('it') # low-voltage activated ca channel
    seg.insert('ihvcn')
    #seg.insert('iH_std')
    seg.insert('leak')
    seg.ena = 10
    seg.ek = -80
    #seg.eh = -40 # Rodrigues and Oertel, 2006

    s = soma()
    if ttx is False:
        s.nav11.gbar = nstomho(1500.0, somaarea) * scalefactor
    else:
        s.nav11.gbar = 0.0 # print s.nav11.gnat
    s.nav11.vsna = 4.3 # was 8
    s.kht.gbar = nstomho(380.0 * 2.0, somaarea) * scalefactor
    s.ka.gbar = nstomho(90.0, somaarea) * scalefactor # was 280
    
    #s.it.gbar = nstomho(14.0 * 4.0, somaarea) * scalefactor
    # this creates a better rebound! with it
    #s.it.vshift = -16
    
    #  was 16.5to allow the tau shift to be about right so it is not so fast.
    #s.iH_std.gbar = nstomho(100.0, somaarea) * scalefactor
    #s.iH_std.vshift = 1.8
    s.ihvcn.gbar = nstomho(220.0, somaarea) * scalefactor
    s.ihvcn.vshift = 16.0
    s.ihvcn.eh = -43
    s.leak.gbar = nstomho(18.0, somaarea) * scalefactor
    s.leak.e = -61
    vm0 = -60.0
    if not runQuiet:
        if message is None:
            print "<< T-stellate: JSR Stellate Type 1 cell model created >>"
        else:
            print message
    return(soma)


def dstellate(debug=False, ttx=False, message=None):
    """ as a type 2-1 from RM03 """
    soma = h.Section() # one compartment of about 29000 um2
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential

    cm = 1.0
    scalefactor = 1.0 # This determines the relative size of the cell
    rinsf = 1.0           # input resistance adjustment (also current...)
    totcap = scalefactor * 25.0
    #  Larger cap value of 25 best fit for mouse D-stellate cell
    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd

    seg = soma
    seg.insert('kht')
    seg.insert('klt')
    seg.insert('nav11')
    seg.insert('ka')
    seg.insert('ihvcn')
    #seg.insert('iH_std')
    seg.insert('leak')
    seg.ena = 20
    seg.ek = -77
    #seg.eh = -43 # Rodrigues and Oertel, 2006
    s = soma()
    if not ttx:
        s.nav11.gbar = nstomho(1000.0, somaarea) * scalefactor
    else:
        s.nav11.gbar = 0.0
    s.nav11.vsna = 8 # voltage shift
    s.kht.gbar = nstomho(250.0, somaarea) * scalefactor
    s.klt.gbar = nstomho(35.0, somaarea) * scalefactor
    s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
    s.ihvcn.gbar = nstomho(3.5, somaarea) * scalefactor
    s.ihvcn.vshift = 0
    s.leak.gbar = nstomho(2, somaarea) * scalefactor
    vm0 = -64.1
    if not runQuiet:
        if message is None:
            print ("<< D-stellate: JSR Stellate",
            " Type I-II cell model created >>")
        else:
            print message
    return(soma)

#
# Integrate and fire version of the d=stellate cell.


def dstellateIF(debug=False, ttx=False, message=None):
    soma = h.Section() # one compartment of about 29000 um2
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential

    cm = 1.0
    scalefactor = 1.0 # This determines the relative size of the cell
    rinsf = 1.0           # input resistance adjustment (also current...)
    totcap = scalefactor * 25.0
    #  Larger value of 25 best fit for mouse D-stellate cell
    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd

    seg = soma
    seg.insert('kht')
    seg.insert('klt')
    seg.insert('nav11')
    seg.insert('ka')
    #seg.insert('ihvcn')
    seg.insert('iH_std')
    seg.insert('leak')
    seg.ena = 20
    seg.ek = -77
    seg.eh = -40 # Rodrigues and Oertel, 2006
    s = soma()
    if not ttx:
        s.nav11.gbar = nstomho(3500.0, somaarea) * scalefactor
    else:
        s.na.gbar = 0.0
    s.nav11.vsna = 8 # voltage shift
    s.kht.gbar = nstomho(150.0, somaarea) * scalefactor
    s.klt.gbar = nstomho(10.0, somaarea) * scalefactor
    s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
    s.iH_std.gbar = nstomho(120.0, somaarea) * scalefactor
    s.iH_std.vshift = 16
    s.lek.gbar = nstomho(24, somaarea) * scalefactor
    vm0 = -64.1
    if not runQuiet:
        if message is None:
            print("<< D-stellate: JSR Stellate",
            " Type I-II cell model created >>")
        else:
            print message
    return(soma)


def dstellate_eager(debug=False, ttx=False, message=None):
    """
    This is a model of the D-Stellate cells as proposed by
    Eager, M.A., Grayden, D.B., Burkitt, A.N., and Meffin, H.,
    "A neural circuit model of the ventral cochlear nucleus",
    Internet:
    http://citeseerx.ist.pus.edu/viewdoc/download?doi=10.1.79.9620.pdf&rep
    =rep&type=pdf
    also cited as:
    Proceedings of theh 10th AUstralian International Conference on
    Speech Science and Technology,
    pp. 539-544, 2004.
    it is based on the Rothman and Manis (2003c) model,
    with small modifications.
    Their model includes dendrites...
    """
    soma = h.Section() # one compartment of about 29000 um2
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential

    cm = 0.9
    Ra = 150
    lstd = 25.0

    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd

    axon = h.Section()
    axon.nseg = 2
    axon.diam = 3.0
    axon.L = 70.0

    axon.connect(soma, 0, 0)

    nDendrite = 2
    dend = nDendrite * [None]

    for i in range(len(dend)):
        dend[i] = h.Section()
        dend[i].nseg = 5
        dend[i].L = 1100.0
        dend[i].diam = 3.5
        dend[i].Ra = 1500.0
        dend[i].cm = 0.9
        dend[i].connect(soma, 0, 1)

    activecompartments = [soma, axon, dend[0], dend[1]]
    #print activecompartments
    for seg in activecompartments:
        s = seg()
        seg.insert('kht')
        seg.insert('klt')
        seg.insert('nav11')
        seg.insert('ka')
        seg.insert('ihvcn')
        #seg.insert('iH_std')
        seg.insert('leak')
        seg.ena = 20
        seg.ek = -77
        # seg.eh = -40 # Rodrigues and Oertel, 2006
        seg.Ra = 150.0
        seg.cm = 0.9
        if not ttx:
            gna = 0.5
        else:
            gna = 0.0
        if seg not in dend:
            s.nav11.gbar = gna # nstomho(gna, somaarea) * scalefactor
            s.nav11.vsna = 8 # voltage shift
            s.kht.gbar = 0.01 # nstomho(500.0, somaarea) * scalefactor
            s.klt.gbar = 0.005 # nstomho(125.0, somaarea) * scalefactor
            s.ka.gbar = 0.0 # nstomho(0.0, somaarea) * scalefactor
            s.ihvcn.gbar = 0.0001 # nstomho(5.0, somaarea) * scalefactor
            s.ihvcn.vshift = 0 # 16
            s.lek.gbar = 0.00025 # nstomho(12.5, somaarea) * scalefactor
        else:
            s.nav11.gbar = 0 # nstomho(gna, somaarea) * scalefactor
            s.nav11.vsna = 0 # voltage shift
            s.kht.gbar = 0.00 # nstomho(500.0, somaarea) * scalefactor
            s.klt.gbar = 0.001 # nstomho(125.0, somaarea) * scalefactor
            s.ka.gbar = 0.0 # nstomho(0.0, somaarea) * scalefactor
            s.ihvcn.gbar = 0.0001 # nstomho(5.0, somaarea) * scalefactor
            s.ihvcn.vshift = 0 # 16
            s.lek.gbar = 0.00025 # nstomho(12.5, somaarea
        print 'nav: ', s.nav11.gbar
        print 'khe: ', s.kht.gbar
        print 'klt: ', s.klt.gbar
        print 'ih:  ', s.ihvcn.gbar
        print 'gleak:', s.leak.g
    vm0 = -64.1
    if not runQuiet:
        if message is None:
            print ("<< D-stellate: Eager et al. ",
            " Type I-II (D-stellate) cell model created >>")
        else:
            print message
    return(soma, dend)


def make_pulse(stim, pulsetype="square"):

    delay = int(numpy.floor(stim['delay'] / h.dt))
    ipi = int(numpy.floor((1000.0 / stim['Sfreq']) / h.dt))
    pdur = int(numpy.floor(stim['dur'] / h.dt))
    posttest = int(numpy.floor(stim['PT'] / h.dt))
    ndur = 5
    if stim['PT'] == 0:
        ndur =1
    maxt = h.dt * (stim['delay'] + (ipi * (stim['NP'] + 3)) +
        posttest + pdur * ndur)
    if 'hold' in stim.keys():
        hold = stim['hold']
    else:
        hold = None
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


def run_iv(ivrange, cell, durs=None, sites=None,
    scales=None, reppulse=None):
    try:
        (imin, imax, istep) = ivrange # unpack the tuple...
    except:
        print ("Cells.py: run_iv argument 1 must have 3 values:",
        " imin, imax and istep")
        sys.exit(1)
    (imin, imax, istep) = ivrange # unpack the tuple...
    #print "min max step: ", imin, imax, istep
    if durs is None:
        durs = [10.0, 100.0, 50.0]
    icur = []
    if reppulse is None:
        istim = h.IClamp2(0.5, sec=cell) # use our new iclamp method
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
        istim = h.iStim(0.5, sec=cell)
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

    #print 'icur: ', icur
    s = cell()

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
        vec['v_soma'].record(cell(0.5)._ref_v, sec=cell)
        vec['ik'].record(cell(0.5)._ref_ik, sec=cell)
        natFlag = False
        try:
            vec['inat'].record(cell(0.5)._ref_inat, sec=cell)
            natFlag = True
        except:
            vec['ina'].record(cell(0.5)._ref_ina, sec=cell)
            pass
        if sites is not None:
            for j in range(len(sites)):
                if sites[j] is not None:
                    print 'section %d : ' % (j),
                    print sites[j]
                    vec['v_meas_%d' % (j)].record(
                        sites[j](0.5)._ref_v, sec=sites[j])
        vec['i_inj'].record(istim._ref_i, sec=cell)
        vec['gh'].record(s.ihvcn._ref_i, sec=cell)
        vec['time'].record(h._ref_t)
        if reppulse is not None:
            vec['i_stim'].play(istim._ref_i, h.dt, 0, sec=cell)

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
    vstim = h.SEClamp(0.5, cell) # use our new iclamp method
    vstim.dur1 = 50.0
    vstim.amp1 = -60
    vstim.dur2 = 500.0
    vstim.amp2 = -60.0
    vstim.dur3 = 400
    vstim.amp3 = -60.0
    vstim.rs = 0.01
    cell.cm = 0.001
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
        vec['v_soma'].record(cell(0.5)._ref_v, sec=cell)
        vec['i_inj'].record(vstim._ref_i, sec=cell)
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



if __name__ == "__main__":
    import argparse
    import sys
    debugFlag = True
    parser = argparse.ArgumentParser(description=('Cells.py:',
    ' Biophysical representatoins of neuorns (mostly auditory)'))
    cclamp = False
    cellinfo = {'types': ['bushy', 'stellate', 'steldend', 'dstellate', 'sgc',
                            'cartwheel', 'pyramidal', 'octopus'],
                'configs': ['std, ''waxon', 'dendrite'],
                'nav': ['std', 'jsrna', 'nav11'],
                'species': ['guineapig', 'guineapig-bushy-II-I',
                                    'guineapig-bushy-II', 'rat', 'mouse'],
                'pulse': ['step', 'pulse']}
    ccivrange = {'bushy': (-1.0, 1.0, 0.1),
                'stellate': (-1.0, 1.0, 0.1),
                'steldend': (-1.0, 1.0, 0.1),
                'dstellate': (-0.25, 1.0, 0.05),
                'sgc:': (-0.5, 0.5, 0.05),
                'cartwheel': (-0.5, 0.5, 0.05),
                'pyramidal': (-1., 1., 0.1),
                'octopus': (-2., 2., 0.2)}
    # scales holds some default scalint to use in the cciv plots
    # argument is {cellname: (xmin, xmax, IVymin, IVymax, FIspikemax,
    # offset(for spikes), crossing (for IV) )}
    ## the "offset" refers to setting the axes back a bit
    scale = {'bushy': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                    'crossing', [0, -60]),
                'stellate': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                    'crossing', [0, -60]),
                'steldend': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                    'crossing', [0, -60]),
                'dstellate': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                    'crossing', [0, -60]),
                'sgc:': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                    'crossing', [0, -60]),
                'cartwheel': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                    'crossing', [0, -60]),
                'pyramidal': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                    'crossing', [0, -60]),
                'octopus': (-1.0, -160., 1.0, -40, 0, 40, 'offset', 5,
                    'crossing', [0, -60])}
    ax = None
    h.celsius = 32
    parser.add_argument('celltype', action='store')
    parser.add_argument('species', action='store', default='guineapig')
     # species is an optional option....
    parser.add_argument('-c', action="store", dest="configuration",
        default='std', help=("Set axon config: %s " %
         [cfg for cfg in cellinfo['configs']]))
    parser.add_argument('--nav', action="store", dest="nav", default="std",
        help=("Choose sodium channel: %s " % [ch for ch in cellinfo['nav']]))
    parser.add_argument('-p', action="store", dest="pulsetype", default="step",
        help=("Set CCIV pulse to step or repeated pulse"))
    clampgroup = parser.add_mutually_exclusive_group()
    clampgroup.add_argument('--vc', action='store_true',
        help="Run in voltage clamp mode")
    clampgroup.add_argument('--cc', action='store_true',
        help="Run in current clamp mode")
    clampgroup.add_argument('--demo', action='store_true',
        help="Run in  voltage clamp demo")
    args = parser.parse_args()
    print args.celltype
    if args.celltype in cellinfo['types']:
        print 'cell: %s is ok' % args.celltype
    else:
        print 'cell: %s is not in our list of cell types' % (args.celltype)
        print 'celltypes: ', cellinfo['types']
        sys.exit(1)

    path = os.path.dirname(__file__)
    h.nrn_load_dll(os.path.join(path, 'i386/special'))
    h.load_file("stdrun.hoc")
    h.load_file(os.path.join(path, "custom_init.hoc"))
    # replace init with one that gets closer to steady state

    print 'configuration: ', args.configuration
    sites = None
    if args.pulsetype == 'step':
        ptype = None
    else:
        ptype = 'pulses'
    if args.configuration in cellinfo['configs']:
        print 'Configuration %s is ok' % args.configuration
    if args.celltype == 'sgc':
        (cell, sgcaxon) = sgc(debugFlag = debugFlag, species = 'mouse',
        nach = 'nav11', chlist = ['ih'])
    elif (args.celltype == 'stellate' and args.nav == 'nav11'
            and args.species == 'guineapig'):
        (cell, dend, axon) = tstellate_rothman_nav11(debug=debugFlag)
    elif (args.celltype == 'steldend'):
        (cell, dendrites, axon) = tstellate_rothman_nav11(debug=debugFlag, dend=True,
         ttx=False, cs=False)
    elif (args.celltype == 'stellate' and args.nav == 'nav11'
            and args.species == 'mouse'):
        cell = tstellate_rothman(species=args.species,
            nav11=True, debug=debugFlag)
    elif (args.celltype == 'bushy' and args.configuration == 'waxon'):
        (cell, [bwaiseg, bwaaxn, bwainternode]) = bushy_waxon(debug=debugFlag)
    elif args.celltype == 'bushy' and args.configuration == 'std':
        (cell, [biseg, baxn, binternode]) = bushy(debug=debugFlag)
        sites = [cell, None, None, None]
    elif args.celltype == 'bushy' and args.configuration == 'dendrite':
        dendriteFlag = True
        (cell, bdends) = bushy(debug=debugFlag, dendrite=dendriteFlag)
        sites = [cell, bdends[0], bdends[1][0]]
    elif args.celltype == 'stellate' and args.nav == 'std':
        cell = tstellate_f(debug=debugFlag)
    elif args.celltype == 'dstellate':
        cell = dstellate(debug=debugFlag)
    else:
        print ("Cell Type %s and configurations nav=%s or config=%s are not available" % (args.celltype, args.nav, args.configuration))
        sys.exit(1)
#    seg = cell()
#
# define the current clamp electrode and default settings
#
    if args.cc is True:
        run_iv(ccivrange[args.celltype], cell,
            scales=scale[args.celltype], sites=sites, reppulse=ptype)
    elif args.vc is True:
        run_vc(-120, -40, 5, cell)
    else:
        if args.demo is True:
            run_democlamp(cell, dendrites)
#-----------------------------------------------------------------------------
#
# If we call this directly, provide a test with the IV function
# - see below to switch cells
#


