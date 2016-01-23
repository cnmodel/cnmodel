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
            raise ValueError ('Pyramidal model %s is unknown', model)


class PyramidalKanold(Pyramidal, Cell):
    """
    DCN pyramidal cell
    Kanold and Manis, 1999, 2001, 2005
    """
    def __init__(self,  morphology=None, decorator=None, hocReader=None, nach='napyr', ttx=False,
                debug=False, species='rat', modelType=None):
        """
        initialize a pyramidal cell, based on the Kanold-Manis (2001) pyramidal cell model.
        Modifications to the cell can be made by calling methods below. These include:
            Converting to a model with modified size and conductances (experimental).
        
        Parameters
        ----------
        morphology : string (default: None)
            a file name to read the cell morphology from. If a valid file is found, a cell is constructed
            as a cable model from the hoc file.
            If None (default), the only a point model is made, exactly according to RM03.
            
        decorator : Python function (default: None)
            decorator is a function that "decorates" the morphology with ion channels according
            to a set of rules.
            If None, a default set of channels aer inserted into the first soma section, and the
            rest of the structure is "bare".
        
        hocReader : Python function (default: None)
            hocReader is the reader that will be used to parse the morphology file, generate
            and connect NEURON sections for the model. The standard hocReader will be the HocReader
            class from neuronvis.
            
        nach : string (default: 'na')
            nach selects the type of sodium channel that will be used in the model. A channel mechanims
            by that name must exist. 
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
        
        species: string (default 'guineapig')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored.
            
        modelType: string (default: None)
            modelType specifies the type of the model that will be used (e.g., "II", "II-I", etc).
            modelType is passed to the decorator, or to species_scaling to adjust point models.
            
        debug: boolean (default: False)
            debug is a boolean flag. When set, there will be multiple printouts of progress and parameters.
            
        Returns
        -------
            Nothing
        
        """
        super(PyramidalKanold, self).__init__()
        if hocReader is not None:
            self.set_reader(hocReader)
        if modelType == None:
            modelType = 'POK'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'Pyramidal',
                       'morphology': morphology, 'decorator': decorator,
                   }

        self.i_test_range=(-0.15, 0.15, 0.01)

        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            soma = h.Section(name="Pyramidal_Soma_%x" % id(self)) # one compartment of about 29000 um2
            soma.nseg = 1
        else:
            """
            instantiate a structured model with the morphology as specified by 
            the morphology file
            """
            soma = self.morphology_from_hoc(morphology=morphology, somasection='sections[0]')

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['napyr', 'kdpyr', 'kif', 'kis', 'ihpyr', 'ihvcn', 'leak', 'kcnq', 'nap']
            for mech in self.mechanisms:
                try:
                    soma.insert(mech)
                except ValueError:
                    print 'WARNING: Mechanism %s not found' % mech
            soma().kif.kif_ivh = -89.6
            self.add_section(soma, 'soma')
            self.species_scaling(silent=True, species=species, modelType=modelType)  # set the default type I-c  cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorated = decorator(self.hr, cellType='Pyramidal', modelType=modelType,
                                 parMap=None)
            self.decorated.channelValidate(self.hr, verify=False)
            self.mechanisms = self.decorated.hf.mechanisms  # copy out all of the mechanisms that were inserted
#        print 'Mechanisms inserted: ', self.mechanisms
        self.get_mechs(soma)
        self.cell_initialize()
        if debug:
            print "<< PYR: POK Pyramidal Cell created >>"

    def species_scaling(self, silent=True, species='rat', modelType='I'):
        soma = self.soma
        if species == 'rat' and modelType in ['I', 'POK']:  # canonical K&M2001 model cell
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

        elif species == 'rat' and modelType == 'II':
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
            raise ValueError('Species %s or species-modelType %s is not recognized for Pyramidal cells' % (species, modelType))

        self.status['species'] = species
        self.status['modelType'] = modelType
        self.cell_initialize(showinfo=True)
        if not silent:
            print 'set cell as: ', species, modelType
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
