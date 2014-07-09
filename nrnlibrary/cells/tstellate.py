from neuron import h
import numpy as np
import neuron as nrn
from ..pynrnutilities import nstomho

from .cell import Cell

__all__ = ['TStellate', 'TStellateNav11', 'TStellateFast'] 

class TStellate(Cell):
    """
    VCN T-stellate model.
    Rothman and Manis, 2003abc (Type I)    
    """
    def __init__(self, nach='nacn', ttx=False, debug=False):
        """
        initialize a planar stellate (T-stellate) cell, using the default parameters for guinea pig from
        R&M2003, as a type I cell.
        Modifications to the cell can be made by calling methods below. These include:
            Converting to a type IA model (add transient K current) (species: guineapig-TypeIA).
            Changing "species" to mouse or cat (scales conductances)
        """
        super(TStellate, self).__init__()

        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': 'guineapig-TypeI', 'ttx': ttx}
        self.e_k = -80  # potassium reversal potential, mV
        self.e_na = 50
        self.c_m = 1.0  # specific membrane capacitance,  uf/cm^2
        self.R_a = 150  # axial resistivity of cytoplasm/axoplasm, ohm.cm
        self.totcap = None
        self.somaarea = None
        self.initsegment = None  # hold initial segment sections
        self.axnode = None  # hold nodes of ranvier sections
        self.internode = None  # hold internode sections
        self.maindend = None  # hold main dendrite sections
        self.secdend = None  # hold secondary dendrite sections
        self.axonsf = None  # axon diameter scale factor
        self.vm0 = -60.

        soma = h.Section() # one compartment of about 29000 um2

        soma.nseg = 1

        if nach == 'nacn':
            soma.insert('nacn')
        elif nach == 'nav11':
            soma.insert('nav11')
        elif nach == 'jsrna':
            soma.insert('jsrna')
        else:
            raise ValueError('Sodium channel %s in type 1 cell not known' % nach)

        soma.insert("kht")
        soma.insert('ka')
        soma.insert('ihvcn')
        soma.insert('leak')
        soma.ek = self.e_k
        self.mechanisms = ['kht', 'ka', 'ihvcn', 'leak', nach]
        self.soma = soma
        self.species_scaling(soma)  # set the default type II cell parameters

        if debug:
                print "<< T-stellate: JSR Stellate Type 1 cell model created >>"
        self.print_mechs(self.soma)

    def set_soma_size_from_Cm(self, cap):
        self.totcap = cap
        self.somaarea = self.totcap * 1E-6 / self.c_m  # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * ((self.somaarea / np.pi) ** 0.5)  # convert from cm to um
        self.soma.diam = lstd
        self.soma.L = lstd

    def species_scaling(self, soma, species='guineapig-TypeI'):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        """
        if species == 'mouse':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            print 'Mouse Tstellate cell'
            self.set_soma_size_from_Cm(25.0)
            self.adjust_na_chans(soma, gbar=800.)
            soma().kht.gbar = nstomho(250.0, self.somaarea)
            soma().ka.gbar = nstomho(0.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(18.0, self.somaarea)
            soma().ihvcn.eh = -43 # Rodrigues and Oertel, 2006
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.vm0 = -60.0
            self.axonsf = 0.5
        if species == 'guineapig-TypeI':  # values from R&M 2003, Type I
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().ka.gbar = nstomho(0.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(0.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.vm0 = -63.6
            self.axonsf = 0.5
        if species == 'guineapig-TypeIt':
            # guinea pig data from Rothman and Manis, 2003, type It
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(80.0, self.somaarea)
            soma().ka.gbar = nstomho(65.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(0.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.vm0 = -64.2
            self.axonsf = 0.5
        if species == 'cat':  # a cat is a big guinea pig Type I
            self.set_soma_size_from_Cm(30.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().ka.gbar = nstomho(0.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(0.5, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.vm0 = -63.6
            self.axonsf = 1.0
        self.status['species'] = species

    def adjust_na_chans(self, soma, gbar=1000., debug=False):
        """
        adjust the sodium channel conductance
        :param soma: a soma object whose sodium channel complement will have it's
        conductances adjusted depending on the channel type
        :return nothing:
        """
        if self.status['ttx']:
            gnabar = 0.0
        else:
            gnabar = nstomho(gbar, self.somaarea)
        nach = self.status['na']
        if nach == 'jsrna':
            soma().jsrna.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'jsrna gbar: ', soma().jsrna.gbar
        elif nach == 'nav11':
            soma().nav11.gbar = gnabar * 0.5
            soma.ena = self.e_na
            soma().nav11.vsna = 4.3
            if debug:
                print "bushy using inva11"
            print 'nav11 gbar: ', soma().nav11.gbar
        else:
            soma().nacn.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'nacn gbar: ', soma().nacn.gbar

    def add_axon(self):
        Cell.add_axon(self, self.soma, self.somaarea, self.c_m, self.R_a, self.axonsf)


class TStellateNav11(Cell):
    """
    VCN T-stellate cell setup from Rothman and Manis, 2003, 
    using nav11 sodium channel model
    
    ttx: if True, turns off sodium channels
    cs: if True, turns off K channels (e.g., cesium in pipette). 
    dend: if True, adds dendrites to the model, based roughly on White et al.,
    1994)
    NOTE: This has been modified from it's original from
    for use in simulating MOUSE stellate cells.
    """
    def __init__(self, debug=False, ttx=False, cs = False, message=None, dend=False):
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
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

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
                ds.leak.gbar = 0.0001
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
        s.leak.gbar = nstomho(2.0, somaarea) * scalefactor
        vm0 = -63.9
        if debug:
            if message is None:
                print ("<< T-stellate: JSR Stellate Type 1 cell model created",
                " - modified for mouse >>")
            else:
                print message
    # print dendrites
        self.soma = soma
        self.dendrites = dendrites

class TStellateFast(Cell):
    """ 
    VCN t-stellate model based on Rothman and Manis 2003, but with fast sodium
    channel 
    """
    def __init__(self, debug=False, ttx=False, message=None, dend=False):
        soma = h.Section() # one compartment of about 29000 um2
        v_potassium = -80       # potassium reversal potential
        v_sodium = 50           # sodium reversal potential

        cm = 1.0
        scalefactor = 1.0 # This determines the relative size of the cell
        rinsf = 1.0           # input resistance adjustment (also current...)
        totcap = scalefactor * 20.0 # cap in pF for cell
        effcap = totcap # sometimes we change capacitance - that's effcap
        somaarea = totcap * 1E-6 / cm # pf -> uF, cm = 1uf/cm^2 nominal;
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

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
        if debug:
            if message is None:
                print "<< T-stellate: JSR Stellate Type 1 cell model created >>"
            else:
                print message
        
        self.soma = soma
