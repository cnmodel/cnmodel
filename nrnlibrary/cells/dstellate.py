import numpy as np
from neuron import h
import neuron as nrn
from ..util import nstomho
from .cell import Cell
from .. import synapses
from .. import data

__all__ = ['DStellate', 'DStellateRothman', 'DStellateEager']


class DStellate(Cell):
    
    type = 'dstellate'

    @classmethod
    def create(cls, model='RM03', **kwds):
        if model == 'RM03':
            return DStellateRothman(**kwds)
        elif model == 'Eager':
            return DStellateEager(**kwds)
        else:
            raise ValueError ('DStellate type %s is unknown', type)

    def make_psd(self, terminal, **kwds):
        from .. import cells

        pre_sec = terminal.section
        pre_cell = terminal.cell
        post_sec = self.soma
        
        if isinstance(pre_cell, cells.SGC):
            # Max conductances for the glu mechanisms are calibrated by 
            # running `synapses/tests/test_psd.py`. The test should fail
            # if these values are incorrect:
            AMPA_gmax = 0.526015135636368
            NMDA_gmax = 0.28738714531937265
            
            # Get AMPAR kinetic constants from database 
            params = data.get('sgc_ampa_kinetics', species='mouse', post_type='dstellate',
                              field=['Ro1', 'Ro2', 'Rc1', 'Rc2', 'PA'])
            
            return synapses.GluPSD(post_sec, terminal,
                                   ampa_gmax=AMPA_gmax,
                                   nmda_gmax=NMDA_gmax,
                                   ampa_params=dict(
                                        Ro1=params['Ro1'],
                                        Ro2=params['Ro2'],
                                        Rc1=params['Rc1'],
                                        Rc2=params['Rc2'],
                                        PA=params['PA']),
                                   **kwds)
        elif isinstance(pre_cell, cells.DStellate):
            # Get GLY kinetic constants from database 
            params = data.get('gly_kinetics', species='mouse', post_type='dstellate',
                              field=['KU', 'KV', 'XMax'])
            psd = synapses.GlyPSD(post_sec, terminal,
                                   psdType='glyfast',
                                   **kwds)
            return psd
        else:
            raise TypeError("Cannot make PSD for %s => %s" % 
                            (pre_cell.__class__.__name__, 
                             self.__class__.__name__))

    def make_terminal(self, post_cell, **kwds):
        from .. import cells
        #
        # set parameters according to the target cell type
        #

        if isinstance(post_cell, cells.Bushy):
            nzones, delay = 10, 0
        elif isinstance(post_cell, cells.TStellate):
            nzones, delay = 5, 0
        elif isinstance(post_cell, cells.DStellate):
            nzones, delay = 5, 0
        else:
            raise NotImplementedError("No knowledge as to how to connect DStellate to cell type %s" %
                                      type(post_cell))
        
        pre_sec = self.soma
        
        return synapses.StochasticTerminal(pre_sec, post_cell, nzones=nzones, 
                                           delay=delay, **kwds)


class DStellateRothman(DStellate):
    """
    VCN D-stellate model:
    as a type I-II from Rothman and Manis, 2003
    """
    def __init__(self, nach='na', ttx=False, debug=False, species='guineapig', type=None):
        """
        initialize a radial stellate (D-stellate) cell, using the default parameters for guinea pig from
        R&M2003, as a type I-II cell.
        Modifications to the cell can be made by calling methods below. These include:
            changing the sodium channel
            Changing "species" to mouse or cat (scales conductances)
            Shifting model type
        """
        super(DStellateRothman, self).__init__()
        print 'rm03 model'
        if type == None:  # allow us to pass None to get the default
            type = 'I-II'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'DStellate'}
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
        soma().leak.erev = self.e_leak
        self.add_section(soma, 'soma')
        self.set_soma_size_from_Cm(12.0)
        self.mechanisms = ['kht', 'klt', 'ihvcn', 'leak', nach]
        self.get_mechs(soma)
        self.species_scaling(species=species, type=type)  # set the default type II cell parameters
        if debug:
                print "<< D-stellate: JSR Stellate Type I-II cell model created >>"

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
        self.status['species'] = species
        self.status['type'] = type
        self.cell_initialize(showinfo=False)
        if not silent:
            print 'set cell as: ', species
            print ' with Vm rest = %6.3f' % self.vm0


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


class DStellateEager(DStellate):
    """
    This is a model of the VCN D-Stellate cells as proposed by
    Eager, M.A., Grayden, D.B., Burkitt, A.N., and Meffin, H.,
    "A neural circuit model of the ventral cochlear nucleus",
    Internet:
    http://citeseerx.ist.pus.edu/viewdoc/download?doi=10.1.79.9620.pdf&rep
    =rep&type=pdf
    also cited as:
    Proceedings of the 10th Australian International Conference on
    Speech Science and Technology, pp. 539-544, 2004.
    It is based on the Rothman and Manis (2003c) model,
    with small modifications.
    Their model includes dendrites and an axon, which are added in this version
    """
    def __init__(self, nach='na', ttx=False, debug=False, species='guineapig', type='I-II'):
        super(DStellateEager, self).__init__()

        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'DStellateEager'}
        self.i_test_range=(-0.25, 0.25, 0.025)  # set range for ic command test

        soma = h.Section(name="DStellateEager_Soma_%x" % id(self)) # one compartment

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
        soma().leak.erev = self.e_leak
        self.mechanisms = ['kht', 'klt', 'ihvcn', 'leak', nach]
        self.add_section(soma, 'soma')
        self.get_mechs(soma)
        self.species_scaling(silent=False, species=species, type=type)  # set the default type II cell parameters
        self.add_axon()  # must follow species scaling so that area parameters are available
        self.add_dendrites()   # similar for dendrites
        if debug:
                print "<< D-stellateEager: Eager DStellate Type I-II cell model created >>"

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
            self.set_soma_size_from_Diam(25.0)
            self.adjust_na_chans(soma, gbar=1000.*0.75)
            soma().kht.gbar = 0.02  # nstomho(150.0, self.somaarea)
            soma().klt.gbar = 0.005  #nstomho(20.0, self.somaarea)
            soma().ihvcn.gbar = 0.0002  #nstomho(2.0, self.somaarea)
            soma().leak.gbar = 0.0005  # nstomho(2.0, self.somaarea)
            self.axonsf = 1.0
        elif species == 'cat' and type == 'I=II':  # a cat is a big guinea pig Type I
            self.set_soma_size_from_Cm(35.0)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(20.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(2.0, self.somaarea)
            soma().leak.gbar = nstomho(2.0, self.somaarea)
            self.axonsf = 1.0
        else:
            raise ValueError('Species %s or species-type %s is not recognized for D-StellateEager cells' %  (species, type))
        self.status['species'] = species
        self.status['type'] = type
        self.cell_initialize(showinfo=True)
        if not silent:
            print ' set cell as: ', species
            print ' with Vm rest = %6.3f' % self.vm0

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
                print 'using jsrna with gbar: ', soma().jsrna.gbar
        elif nach == 'nav11':
            soma().nav11.gbar = gnabar * 0.5
            soma.ena = self.e_na
            soma().nav11.vsna = 4.3
            if debug:
                print "using inva11 with gbar:", soma().na.gbar
            print 'nav11 gbar: ', soma().nav11.gbar
        elif nach == 'na':
            soma().na.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'using na with gbar: ', soma().na.gbar
        elif nach == 'nach':
            soma().nach.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'uwing nacn with gbar: ', soma().nacn.gbar
        else:
            raise ValueError("DstellateEager setting Na channels: channel %s not known" % nach)
        print soma().na.gbar

    def add_axon(self):
        #Cell.add_axon(self, nodes=1, c_m=self.c_m, R_a=self.R_a, axonsf=self.axonsf, dia=3.0, len=70, seg=2)
        # The Eager et al model just uses one cable, 70 microns long and 3 microns in dameter.
        naxons = 1
        axon = []
        for i in range(naxons):
            axon.append(h.Section(cell=self.soma))
        for i in range(naxons):
            axon[i].connect(self.soma)
            axon[i].L = 70
            axon[i].diam = 3.0
            axon[i].Ra = 500
            axon[i].cm = 0.9
            axon[i].nseg = 2
            axon[i].insert("kht")
            axon[i].insert('klt')
            axon[i].insert('ihvcn')
            axon[i].insert('leak')
            axon[i].insert('na')
            axon[i].ek = self.e_k
            axon[i].ena = self.e_na
            axon[i]().leak.erev = self.e_leak
            axon[i]().na.gbar = 0.5
            axon[i]().klt.gbar = 0.005
            axon[i]().kht.gbar = 0.02
            axon[i]().ihvcn.gbar = 0.0002
            axon[i]().leak.gbar = 0.0005
        self.status['axon'] = True
        self.add_section(axon, 'axon')

    def add_dendrites(self):
        """
        The Eager model uses simple passive dendrites, which are built here
        """
        nDend = range(2) # these will be simple, unbranced, N=4 dendrites
        dendrites=[]
        for i in nDend:
            dendrites.append(h.Section(cell=self.soma))
        for i in nDend:
            dendrites[i].connect(self.soma)
            dendrites[i].L = 1100 # length of the dendrite (not tapered)
            dendrites[i].diam = 3.5 # dendrite diameter
            dendrites[i].nseg = 5 # # segments in dendrites
            dendrites[i].Ra = 1500 # ohm.cm
            dendrites[i].insert('leak') # leak
            dendrites[i]().leak.gbar = 0.00025
            dendrites[i]().leak.erev = self.e_leak
        self.maindend = dendrites
        self.status['dendrites'] = True
        self.add_section(self.maindend, 'maindend')


