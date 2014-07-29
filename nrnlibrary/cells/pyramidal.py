from neuron import h
import neuron as nrn
from ..util import nstomho
import numpy as np

from .cell import Cell

__all__ = ['Pyramidal']


class Pyramidal(Cell):
    """
    DCN pyramidal cell
    Kanold and Manis, 1999, 2001, 2005
    """
    def __init__(self, nach='napyr', ttx=False, debug=False, species='rat', type=None):
        """
        initialize a planar stellate (T-stellate) cell, using the default parameters for guinea pig from
        R&M2003, as a type I cell.
        Modifications to the cell can be made by calling methods below. These include:
            Converting to a type IA model (add transient K current) (species: guineapig-TypeIA).
            Changing "species" to mouse or cat (scales conductances)
        """
        super(Pyramidal, self).__init__()
        if type == None:
            type = 'I'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'Pyramidal'}

        self.i_test_range=(-0.15, 0.15, 0.01)

        soma = h.Section(name="Pyramidal_Soma_%x" % id(self)) # one compartment of about 29000 um2

        soma.nseg = 1

        # if nach in ['nacn', 'na']:
        #     soma.insert('na')
        # elif nach == 'nav11':
        #     soma.insert('nav11')
        # elif nach == 'jsrna':
        #     soma.insert('jsrna')
        # else:
        #     raise ValueError('Sodium channel %s in type 1 cell not known' % nach)

        self.mechanisms = ['napyr', 'kdpyr', 'kif', 'kis', 'ihpyr', 'leak']
        for mech in self.mechanisms:
            soma.insert(mech)
        soma().kif.kif_ivh = -89.6
        self.add_section(soma, 'soma')
        self.species_scaling(silent=True, species=species, type=type)  # set the default type I-c  cell parameters
        self.get_mechs(soma)
#        self.cell_initialize()
        if debug:
            print "<< PYR: POK Pyramidal Cell created >>"

    def species_scaling(self, silent=True, species='rat', type='I'):
        soma = self.soma
        if species == 'rat' and type == 'I':
            self.set_soma_size_from_Cm(12.0)
            soma().napyr.gbar = nstomho(350, self.somaarea)
            #soma().pyr.gnapbar = 0.0
            soma().kdpyr.gbar = nstomho(80, self.somaarea) # used to be 20?
            soma().kif.gbar = nstomho(150, self.somaarea)
            soma().kis.gbar = nstomho(40, self.somaarea)
            soma().ihpyr.gbar = nstomho(3, self.somaarea)
            soma().leak.gbar = nstomho(2.8, self.somaarea)
            soma().leak.erev = -57.7  # override default values in cell.py
            soma().ena = 50.0
            soma().ek = -81.5
            soma().eh = -43

        else:
            raise ValueError('Species %s or species-type %s is not recognized for T-stellate cells' % (species, type))

        self.status['species'] = species
        self.status['type'] = type
        self.cell_initialize(showinfo=True)
        if not silent:
            print 'set cell as: ', species
            print ' with Vm rest = %f' % self.vm0

    # more complex cell type:
    #        gkpksk = nstomho(0, self.somaarea)
    #        gkir = nstomho(0, self.somaarea) # incude KIR here, but set to 0

    # # set up soma like a pyramidal cell
    #     soma.nseg = 1
    #     soma.diam = lstd
    #     soma.L = lstd # these are expressed in microns...
    #     soma.insert('pyr')
    #     soma.insert('kpksk')
    #     soma.insert('cadiff') # diffusion
    #     soma.insert('cap') # p-type calcium current
    #     #soma.insert('nacum') # sodium accumulation (yes!)
    #     soma.insert('nakpump') # and of course a pump to handle it.
    #     soma.insert('k_conc')
    #     soma.insert('na_conc')
    #     # soma.insert('kna')
    #
    #     seg = soma()
    #     seg.kpksk.gbar = gkpksk
    #     seg.cap.pcabar = 0.00002
    #     seg.pyr.gbar = gnab
    #     seg.pyr.gbar = gnap
    #     seg.pyr.gbar = gkb
    #     seg.pyr.gbar = gkfb
    #     seg.pyr.gbar = gksb
    #     seg.pyr.gl = glb
    #     seg.pyr.gbar = ghb
    # # seg.pyr.gbar = gkir
    #     seg.pyr.kif_ivh = -89.6
    #

    def i_currents(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point
        vrange brackets the interval
        Overrides i_currents in cells.py
        """
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = V
        #h.finitialize()
        self.ix = {}

        if 'kif' in self.mechanisms:
             self.ix['kif'] = self.soma().kif.gkif*(V - self.soma().ek)
        if 'kis' in self.mechanisms:
             self.ix['kis'] = self.soma().kis.gkis*(V - self.soma().ek)
        if 'ihpyr' in self.mechanisms:
             self.ix['ihpyr'] = self.soma().ihpyr.gh*(V - self.soma().eh)
        if 'napyr' in self.mechanisms:
             self.ix['napyr'] = self.soma().napyr.gna*(V - self.soma().ena)
        if 'kdpyr' in self.mechanisms:
             self.ix['kdpyr'] = self.soma().kdpyr.gk*(V - self.soma().ek)

        # leak
        if 'leak' in self.mechanisms:
            self.ix['leak'] = self.soma().leak.gbar*(V - self.soma().leak.erev)
        return np.sum([self.ix[i] for i in self.ix])
