from neuron import h
import neuron as nrn
from ..util import nstomho
import numpy as np

from .cell import Cell

__all__ = ['Pyramidal', 'PyramidalKanold']

class Pyramidal(Cell):

    @classmethod
    def create(cls, model='POK', **kwds):
        if model == 'POK':
            return PyramidalKanold(**kwds)
        else:
            raise ValueError ('DStellate type %s is unknown', type)


class PyramidalKanold(Pyramidal, Cell):
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
        super(PyramidalKanold, self).__init__()
        if type == None:
            type = 'I'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'type': type, 'ttx': ttx, 'name': 'Pyramidal'}

        self.i_test_range=(-0.15, 0.15, 0.01)

        soma = h.Section(name="Pyramidal_Soma_%x" % id(self)) # one compartment of about 29000 um2

        soma.nseg = 1

        self.mechanisms = ['napyr', 'kdpyr', 'kif', 'kis', 'ihpyr', 'ihvcn', 'leak', 'kcnq', 'nap']
        for mech in self.mechanisms:
            try:
                soma.insert(mech)
            except ValueError:
                print 'WARNING: Mechanism %s not found' % mech
        soma().kif.kif_ivh = -89.6
        self.add_section(soma, 'soma')
        self.species_scaling(silent=False, species=species, type=type)  # set the default type I-c  cell parameters
        self.get_mechs(soma)
#        self.cell_initialize()
        if debug:
            print "<< PYR: POK Pyramidal Cell created >>"

    def species_scaling(self, silent=True, species='rat', type='I'):
        soma = self.soma
        if species == 'rat' and type == 'I':  # canonical K&M2001 model cell
            self.set_soma_size_from_Cm(12.0)
            soma().napyr.gbar = nstomho(350, self.somaarea)
            soma().nap.gbar = 0.0
            soma().kdpyr.gbar = nstomho(80, self.somaarea) # used to be 20?
            soma().kcnq.gbar = 0 # does not exist in canonical model.
            soma().kif.gbar = nstomho(150, self.somaarea)
            soma().kis.gbar = nstomho(40, self.somaarea)
            soma().ihpyr.gbar = nstomho(2.8, self.somaarea)
            soma().ihvcn.gbar = nstomho(0., self.somaarea)
            soma().leak.gbar = nstomho(2.8, self.somaarea)
            soma().leak.erev = -62  # override default values in cell.py
            soma().ena = 50.0
            soma().ek = -81.5
            soma().ihpyr.eh = -43
            soma().ihvcn.eh = -43

        elif species == 'rat' and type == 'II':
            """
            Modified canonical K&M2001 model cell
            In this model version, the specific membrane capacitance is modified
            so that the overall membrane time constant is consistent with experimental
            measures in slices. However, this is not a physiological value. Attempts
            to use the normal 1 uF/cm2 value were unsuccessful in establishing the expected
            ~12 msec time constant.
            This model also adds a KCNQ channel, as described by Li et al., 2012.
            """
            self.c_m = 20.0
            self.set_soma_size_from_Diam(30.0)
            #self.set_soma_size_from_Cm(80.0)
            print 'diameter: %7.1f' % self.soma.diam
            self.refarea = self.somaarea
            soma().napyr.gbar = nstomho(550, self.refarea)
            soma().nap.gbar = nstomho(20.0, self.refarea)
            soma().kcnq.gbar = nstomho(2, self.refarea)  # pyramidal cells have kcnq: Li et al, 2011 (Thanos)
            soma().kdpyr.gbar = nstomho(180, self.refarea) # Normally 80.
            soma().kif.gbar = nstomho(150, self.refarea) # normally 150
            soma().kis.gbar = nstomho(120, self.refarea) # 40
            soma().ihpyr.gbar = nstomho(2.8, self.refarea)
            soma().ihvcn.gbar = nstomho(0., self.refarea)
            soma().leak.gbar = nstomho(1.5, self.refarea)
            soma().leak.erev = -62.  # override default values in cell.py
            soma().ena = 50.0
            soma().ek = -81.5
            soma().ihpyr.eh = -43
            soma().ihvcn.eh = -43
            if not self.status['dendrites']:
                self.add_dendrites()

        else:
            raise ValueError('Species %s or species-type %s is not recognized for Pyramidal cells' % (species, type))

        self.status['species'] = species
        self.status['type'] = type
        self.cell_initialize(showinfo=True)
        if not silent:
            print 'set cell as: ', species, type
            print ' with Vm rest = %f' % self.vm0
            print self.status
            for m in self.mechanisms:
                print '%s.gbar = %f' % (m, eval('soma().%s.gbar' % m))

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
        h.finitialize()
        self.ix = {}

        if 'napyr' in self.mechanisms:
             self.ix['napyr'] = self.soma().napyr.gna*(V - self.soma().ena)
        if 'nap' in self.mechanisms:
             self.ix['nap'] = self.soma().nap.gnap*(V - self.soma().ena)
        if 'kdpyr' in self.mechanisms:
             self.ix['kdpyr'] = self.soma().kdpyr.gk*(V - self.soma().ek)
        if 'kif' in self.mechanisms:
             self.ix['kif'] = self.soma().kif.gkif*(V - self.soma().ek)
        if 'kis' in self.mechanisms:
             self.ix['kis'] = self.soma().kis.gkis*(V - self.soma().ek)
        if 'kcnq' in self.mechanisms:
             self.ix['kcnq'] = self.soma().kcnq.gk*(V - self.soma().ek)
        if 'ihvcn' in self.mechanisms:
             self.ix['ihvcn'] = self.soma().ihvcn.gh*(V - self.soma().ihvcn.eh)
        if 'ihpyr' in self.mechanisms:
             self.ix['ihpyr'] = self.soma().ihpyr.gh*(V - self.soma().ihpyr.eh)
        # leak
        if 'leak' in self.mechanisms:
            self.ix['leak'] = self.soma().leak.gbar*(V - self.soma().leak.erev)
        return np.sum([self.ix[i] for i in self.ix])

    def add_dendrites(self):
        """
        Add simple unbranched dendrite.
        The dendrites have some kd, kif and ih current
        """
        nDend = range(2) # these will be simple, unbranced, N=4 dendrites
        dendrites=[]
        for i in nDend:
            dendrites.append(h.Section(cell=self.soma))
        for i in nDend:
            dendrites[i].connect(self.soma)
            dendrites[i].L = 250 # length of the dendrite (not tapered)
            dendrites[i].diam = 1
            dendrites[i].cm = self.c_m
            #h('dendrites[i].diam(0:1) = 2:1') # dendrite diameter, with tapering
            dendrites[i].nseg = 21 # # segments in dendrites
            dendrites[i].Ra = 150 # ohm.cm
            dendrites[i].insert('napyr')
            dendrites[i]().napyr.gbar = 0.00
            dendrites[i].insert('kdpyr')
            dendrites[i]().kdpyr.gbar = 0.002 # a little Ht
            dendrites[i].insert('kif')
            dendrites[i]().kif.gbar = 0.0001 # a little Ht
            dendrites[i].insert('leak') # leak
            dendrites[i]().leak.gbar = 0.00001
            dendrites[i].insert('ihvcn') # some H current
            dendrites[i]().ihvcn.gbar =  0. # 0.00002
            dendrites[i]().ihvcn.eh = -43.0
        self.maindend = dendrites
        self.status['dendrites'] = True
        self.add_section(self.maindend, 'maindend')
