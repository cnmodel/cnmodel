from neuron import h
import neuron as nrn
from ..util import nstomho
import numpy as np
from .cell import Cell

__all__ = ['DStellate', 'DStellateIF', 'DStellateEager']


class DStellate(Cell):
    """ 
    VCN D-stellate model:
    as a type I-II from Rothman and Manis, 2003
    """
    def __init__(self, nach='na', ttx=False, debug=False):
        """
        initialize a radial stellate (D-stellate) cell, using the default parameters for guinea pig from
        R&M2003, as a type I-II cell.
        Modifications to the cell can be made by calling methods below. These include:
            changing the sodium channel
            Changing "species" to mouse or cat (scales conductances)
            Shifting model type
        """
        super(DStellate, self).__init__()

        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': 'guineapig', 'type': 'I-II', 'ttx': ttx, 'name': 'DStellate'}
        self.e_k = -70  # potassium reversal potential, mV
        self.e_na = 50
        self.c_m = 0.9  # specific membrane capacitance,  uf/cm^2
        self.R_a = 150  # axial resistivity of cytoplasm/axoplasm, ohm.cm
        self.vm0 = -64.1
        self.i_test_range=(-0.25, 0.25, 0.025)  # set range for ic command test

        soma = h.Section(name="DStellate_Soma_%x" % id(self)) # one compartment

        soma.nseg = 1

        if nach in ['nacn', 'na']:
            soma.insert('na')
        elif nach == 'nav11':
            soma.insert('nav11')
        elif nach == 'jsrna':
            soma.insert('jsrna')
        else:
            raise ValueError('Sodium channel %s in type 1 cell not known' % nach)

        soma.insert("kht")
        soma.insert('klt')
        soma.insert('ihvcn')
        soma.insert('leak')
        soma.ek = self.e_k
        soma().leak.erev = -65
        self.mechanisms = ['kht', 'klt', 'ihvcn', 'leak', nach]
        self.add_section(soma, 'soma')
        self.get_mechs(soma)
        self.cell_initialize(showinfo=True)
        self.species_scaling(silent=False)  # set the default type II cell parameters
        if debug:
                print "<< D-stellate: JSR Stellate Type I-II cell model created >>"

    def set_soma_size_from_Cm(self, cap):
        self.totcap = cap
        self.somaarea = self.totcap * 1E-6 / self.c_m  # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * ((self.somaarea / np.pi) ** 0.5)  # convert from cm to um
        self.soma.diam = lstd
        self.soma.L = lstd

    def species_scaling(self, species='guineapig', type='I-II', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        """
        soma = self.soma
        if species == 'mouse' and type == 'I-II':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            self.set_soma_size_from_Cm(25.0)
            self.adjust_na_chans(soma, gbar=800.)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(20.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(2.0, self.somaarea)
            soma().ihvcn.eh = -43 # Rodrigues and Oertel, 2006
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.axonsf = 0.5
        elif species == 'guineapig' and type == 'I-II':  # values from R&M 2003, Type II-I
            self.set_soma_size_from_Cm(12.0)
            self.adjust_na_chans(soma, gbar=1000.)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(20.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(2.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.axonsf = 0.5
        elif species == 'cat' and type == 'I=II':  # a cat is a big guinea pig Type I
            self.set_soma_size_from_Cm(35.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(20.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(2.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.axonsf = 1.0
        else:
            raise ValueError('Species %s or species-type %s is not recognized for D-Stellate cells' %  (species, type))
        #self.cell_initialize(showinfo=True)
        self.vm0 = self.find_i0()
        if not silent:
            print 'set cell as: ', species
            print ' with Vm rest = %6.3f' % self.vm0

        self.status['species'] = species
        self.status['type'] = type

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
        elif nach == 'na':
            soma().na.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'na gbar: ', soma().na.gbar
        elif nach == 'nach':
            soma().nach.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'nacn gbar: ', soma().nacn.gbar
        else:
            raise ValueError("Dstellate setting Na channels: channel %s not known" % nach)

    def add_axon(self):
        Cell.add_axon(self, self.soma, self.somaarea, self.c_m, self.R_a, self.axonsf)

    def add_dendrites(self):
        """
        Add simple unbranched dendrites to basic Rothman Type I-II model.
        The dendrites have some kht and ih current
        """
        cs = False  # not implemented outside here - internal Cesium.
        nDend = range(4) # these will be simple, unbranced, N=4 dendrites
        dendrites=[]
        for i in nDend:
            dendrites.append(h.Section(cell=self.soma))
        for i in nDend:
            dendrites[i].connect(self.soma)
            dendrites[i].L = 300 # length of the dendrite (not tapered)
            dendrites[i].diam = 1.25 # dendrite diameter
            dendrites[i].nseg = 21 # # segments in dendrites
            dendrites[i].Ra = 150 # ohm.cm
            dendrites[i].insert('kht')
            if cs is False:
                dendrites[i]().kht.gbar = 0.005 # a little Ht
            else:
                dendrites[i]().kht.gbar = 0.0
            dendrites[i].insert('leak') # leak
            dendrites[i]().leak.gbar = 0.0001
            dendrites[i].insert('ihvcn') # some H current
            dendrites[i]().ihvcn.gbar = 0.# 0.001
            dendrites[i]().ihvcn.eh = -43.0
        self.maindend = dendrites
        self.status['dendrites'] = True
        self.add_section(self.maindend, 'maindend')


class DStellateIF(Cell):
    """
    Integrate and fire version of the VCN D-stellate cell.
    """
    def __init__(self, debug=False, ttx=False, message=None):
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
        lstd = 1E4 * ((somaarea / 3.14159) ** 0.5) # convert from cm to um

        soma.nseg = 1
        soma.diam = lstd
        soma.L = lstd

        seg = soma
        seg.insert('kht')
        seg.insert('klt')
        seg.insert('nav11')
        #seg.insert('ka')
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
        #s.ka.gbar = nstomho(0.0, somaarea) * scalefactor
        s.iH_std.gbar = nstomho(120.0, somaarea) * scalefactor
        s.iH_std.vshift = 16
        s.leak.gbar = nstomho(24, somaarea) * scalefactor
        vm0 = -64.1
        if debug:
            if message is None:
                print("<< D-stellate: JSR Stellate",
                " Type I-II cell model created >>")
            else:
                print message
                
        self.soma = soma


class DStellateEager(Cell):
    """
    This is a model of the VCN D-Stellate cells as proposed by
    Eager, M.A., Grayden, D.B., Burkitt, A.N., and Meffin, H.,
    "A neural circuit model of the ventral cochlear nucleus",
    Internet:
    http://citeseerx.ist.pus.edu/viewdoc/download?doi=10.1.79.9620.pdf&rep
    =rep&type=pdf
    also cited as:
    Proceedings of the 10th Australian International Conference on
    Speech Science and Technology,
    pp. 539-544, 2004.
    It is based on the Rothman and Manis (2003c) model,
    with small modifications.
    Their model includes dendrites and an axon, which are added in this version
    """
    def __init__(self, debug=False, ttx=False, message=None):
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
            #seg.insert('ka')
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
                #s.ka.gbar = 0.0 # nstomho(0.0, somaarea) * scalefactor
                s.ihvcn.gbar = 0.0001 # nstomho(5.0, somaarea) * scalefactor
                s.ihvcn.vshift = 0 # 16
                s.leak.gbar = 0.00025 # nstomho(12.5, somaarea) * scalefactor
            else:
                s.nav11.gbar = 0 # nstomho(gna, somaarea) * scalefactor
                s.nav11.vsna = 0 # voltage shift
                s.kht.gbar = 0.00 # nstomho(500.0, somaarea) * scalefactor
                s.klt.gbar = 0.001 # nstomho(125.0, somaarea) * scalefactor
                s.ka.gbar = 0.0 # nstomho(0.0, somaarea) * scalefactor
                s.ihvcn.gbar = 0.0001 # nstomho(5.0, somaarea) * scalefactor
                s.ihvcn.vshift = 0 # 16
                s.leak.gbar = 0.00025 # nstomho(12.5, somaarea
            print 'nav: ', s.nav11.gbar
            print 'khe: ', s.kht.gbar
            print 'klt: ', s.klt.gbar
            print 'ih:  ', s.ihvcn.gbar
            print 'gleak:', s.leak.g
        vm0 = -64.1
        if debug:
            if message is None:
                print ("<< D-stellate: Eager et al. ",
                " Type I-II (D-stellate) cell model created >>")
            else:
                print message
                
        self.soma = soma
        self.dend = dend
