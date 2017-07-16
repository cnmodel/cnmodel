from neuron import h
import numpy as np
#import neuron as nrn

from .cell import Cell
from .. import synapses
from ..util import nstomho
from ..util import Params
from .. import data

__all__ = ['Tuberculoventral'] 


class Tuberculoventral(Cell):
    
    type = 'tuberculoventral'

    @classmethod
    def create(cls, modelType='TVmouse', **kwds):
        if modelType in ['TVmouse', 'I']:
            return Tuberculoventral(**kwds)
        else:
            raise ValueError ('Tuberculoventral type %s is unknown', modelType)

    def make_psd(self, terminal, psd_type, **kwds):
        """
        Connect a presynaptic terminal to one post section at the specified location, with the fraction
        of the "standard" conductance determined by gbar.
        The default condition is to try to pass the default unit test (loc=0.5)
        
        Parameters
        ----------
        terminal : Presynaptic terminal (NEURON object)
        
        psd_type : either simple or multisite PSD for bushy cell
        
        kwds: dict of options. Two are currently handled:
        postsize : expect a list consisting of [sectionno, location (float)]
        AMPAScale : float to scale the ampa currents
        
        """
        if 'postsite' in kwds:  # use a defined location instead of the default (soma(0.5)
            postsite = kwds['postsite']
            loc = postsite[1]  # where on the section?
            uname = 'sections[%d]' % postsite[0]  # make a name to look up the neuron section object
            post_sec = self.hr.get_section(uname)  # Tell us where to put the synapse.
        else:
            loc = 0.5
            post_sec = self.soma
        
        if psd_type == 'simple':
            return self.make_exp2_psd(post_sec, terminal, loc=loc)
        elif psd_type == 'multisite':
            if terminal.cell.type == 'sgc':
                # Max conductances for the glu mechanisms are calibrated by 
                # running `synapses/tests/test_psd.py`. The test should fail
                # if these values are incorrect
                self.AMPAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.type, field='AMPAR_gmax')*1e3
                self.NMDAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.type, field='NMDAR_gmax')*1e3
                self.Pr = data.get('sgc_synapse', species=self.species,
                        post_type=self.type, field='Pr')
                # adjust gmax to correct for initial Pr
                self.AMPAR_gmax = self.AMPAR_gmax/self.Pr
                self.NMDAR_gmax = self.NMDAR_gmax/self.Pr
                if 'AMPAScale' in kwds:
                    self.AMPA_gmax = self.AMPA_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    self.NMDA_gmax = self.NMDA_gmax*kwds['NMDAScale']
                return self.make_glu_psd(post_sec, terminal, self.AMPAR_gmax, self.NMDAR_gmax, loc=loc)
            elif terminal.cell.type == 'dstellate':  # WBI input -Voigt, Nelken, Young
                return self.make_gly_psd(post_sec, terminal, type='glyfast', loc=loc)
            elif terminal.cell.type == 'tuberculoventral':  # TV cells talk to each other-Kuo et al.
                return self.make_gly_psd(post_sec, terminal, type='glyfast', loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.type, self.type))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)

    def make_terminal(self, post_cell, term_type, **kwds):
        pre_sec = self.soma
        if term_type == 'simple':
            return synapses.SimpleTerminal(pre_sec, post_cell, 
                                           spike_source=self.spike_source, **kwds)
        elif term_type == 'multisite':
            if post_cell.type == 'bushy':
                nzones, delay = 10, 0
            elif post_cell.type == 'tstellate':
                nzones, delay = 5, 0
            elif post_cell.type == 'tuberculoventral':
                nzones, delay = 2, 0
            elif post_cell.type == 'pyramidal':
                nzones, delay = 5, 0
            else:
                raise NotImplementedError("No knowledge as to how to connect DStellate to cell type %s" %
                                        type(post_cell))
            
            pre_sec = self.soma
            return synapses.StochasticTerminal(pre_sec, post_cell, nzones=nzones, 
                                            delay=delay, **kwds)
        else:
            raise ValueError("Unsupported terminal type %s" % term_type)


class Tuberculoventral(Tuberculoventral):
    """
    Tuberculoventral Neuron (DCN) base model
    Adapted from T-stellate model, using target parameters from Kuo et al. J. Neurophys. 2012
    """
    def __init__(self, morphology=None, decorator=None, nach=None, ttx=False,
                species='mouse', modelType=None, debug=False):
        """
        Initialize a DCN Tuberculoventral cell, using the default parameters for guinea pig from
        R&M2003, as a type I cell.
        Modifications to the cell can be made by calling methods below. These include:
        Converting to a type IA model (add transient K current) (species: guineapig-TypeIA).
        Changing "species" to mouse or cat (scales conductances)
        
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
        super(Tuberculoventral, self).__init__()
        if modelType == None:
            modelType = 'TVmouse'
        if nach == None:
            nach = 'nacncoop'
            
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'Tuberculoventral',
                       'morphology': morphology, 'decorator': decorator, 'temperature': None}

        self.i_test_range = {'pulse': (-0.4, 0.6, 0.02)}
        self.vrange = [-80., -60.]  # set a default vrange for searching for rmp
        
        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            print "<< Tuberculoventral model: Creating point cell >>"
            soma = h.Section(name="Tuberculoventral_Soma_%x" % id(self))  # one compartment of about 29000 um2
            soma.nseg = 1
            self.add_section(soma, 'soma')
        else:
            """
            instantiate a structured model with the morphology as specified by 
            the morphology file
            """
            print "<< Tuberculoventral model: Creating structured cell >>"
            self.set_morphology(morphology_file=morphology)

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['kht', 'ka', 'ihvcn', 'leak', nach]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ek = self.e_k
            self.soma.ena = self.e_na
            self.soma().ihvcn.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.species_scaling(silent=True, species=species, modelType=modelType)  # adjust the default parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
#        print 'Mechanisms inserted: ', self.mechanisms
        
        self.get_mechs(self.soma)
#        self.cell_initialize(vrange=self.vrange)
        if debug:
                print "<< Tuberculoventral cell model created >>"

    def species_scaling(self, species='guineapig', modelType='TVmouse', silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
        
        Parameters
        ----------
        species : string (default: 'guineapig')
            name of the species to use for scaling the conductances in the base point model
            Must be one of mouse, cat, guineapig
        
        modelType: string (default: 'I-c')
            definition of model type from RM03 models, type I-c or type I-t
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        """
        soma = self.soma
        if species == 'mouse' and modelType == 'TVmouse':
            """#From Kuo 150 Mohm, 10 msec tau
            Firing at 600 pA about 400 Hz
            These values from brute_force runs, getting 380 Hz at 600 pA at 35C
            Input resistance and vm is ok, time constnat is short
                *** Rin:       168  tau:       7.8   v:  -68.4
            Attempts to get longer time constant - cannot keep rate up.
            """
            # Adapted from TStellate model type I-c'
            self.vrange=[-80., -58.]
            self.set_soma_size_from_Cm(35.0)
            self._valid_temperatures = (34.,)
            if self.status['temperature'] is None:
                self.set_temperature(34.)
            self.adjust_na_chans(soma, gbar=5800.)
            soma().kht.gbar = nstomho(400.0, self.somaarea) # was 2000
            soma().ka.gbar = nstomho(65.0, self.somaarea)
            soma().ihvcn.gbar = nstomho(2.5, self.somaarea)  # 1.25
            soma().ihvcn.eh = -43 # Rodrigues and Oertel, 2006
            soma().leak.gbar = nstomho(4.5, self.somaarea)  # 5.5
            soma().leak.erev = -72.0
            self.axonsf = 0.5
        else:
            raise ValueError('Species %s or species-type %s is not recognized for Tuberculoventralcells' % (species, type))

        self.status['species'] = species
        self.status['modelType'] = modelType
        self.check_temperature()


    def channel_manager(self, modelType='RM03'):
        """
        This routine defines channel density maps and distance map patterns
        for each type of compartment in the cell. The maps
        are used by the ChannelDecorator class (specifically, it's private
        _biophys function) to decorate the cell membrane.
        
        Parameters
        ----------
        modelType : string (default: 'RM03')
            A string that defines the type of the model. Currently, 3 types are implemented:
            RM03: Rothman and Manis, 2003 somatic densities for guinea pig
            XM13: Xie and Manis, 2013, somatic densities for mouse
            XM13PasDend: XM13, but with only passive dendrites, no channels.
        
        Returns
        -------
        Nothing
        
        Notes
        -----
        
        This routine defines the following variables for the class:
            
            - conductances (gBar)
            - a channelMap (dictonary of channel densities in defined anatomical compartments)
            - a current injection range for IV's (when testing)
            - a distance map, which defines how selected conductances in selected compartments
                will change with distance. This includes both linear and exponential gradients,
                the minimum conductance at the end of the gradient, and the space constant or
                slope for the gradient.
        
        """
        if modelType == 'TVmouse':
            print 'decorate as tvmouse'
#            totcap = 95.0E-12  # Tuberculoventral cell (type I), based on stellate, adjusted for Kuo et al. TV firing
            self.set_soma_size_from_Section(self.soma)
            totcap = self.totcap
            refarea = self.somaarea # totcap / self.c_m  # see above for units
            self.gBar = Params(nabar=1520.0E-9/refarea,
                               khtbar=160.0E-9/refarea,
                               kltbar=0.0E-9/refarea,
                               kabar=65.0/refarea,
                               ihbar=1.25E-9/refarea,
                               leakbar=5.5E-9/refarea,
            )
            self.channelMap = {
                'axon': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar,
                         'ihvcn': 0., 'leak': self.gBar.leakbar / 4.},
                'hillock': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
                             'ihvcn': 0., 'leak': self.gBar.leakbar, },
                'initseg': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
                            'ihvcn': self.gBar.ihbar / 2.,
                            'leak': self.gBar.leakbar, },
                'soma': {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar,
                         'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar,
                         'leak': self.gBar.leakbar, },
                'dend': {'nacn': self.gBar.nabar / 2.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.5,
                         'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
                'apic': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.2,
                         'ihvcn': self.gBar.ihbar / 4.,
                         'leak': self.gBar.leakbar * 0.2, },
            }
            self.irange = np.linspace(-0.3, 0.6, 10)
            self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                            'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                     'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # gradients are: flat, linear, exponential
                            }
        else:
            raise ValueError('model type %s is not implemented' % modelType)
        
    def adjust_na_chans(self, soma, gbar=1000., debug=False):
        """
        Adjust the sodium channel conductance, depending on the type of conductance
        
        Parameters
        ----------
        soma : NEURON section object (required)
            This identifies the soma object whose sodium channel complement will have it's
            conductances adjusted depending on the sodium channel type
        gbar : float (default: 1000.)
            The "maximal" conductance to be set in the model.
        debug : boolean (default: False)
            A flag the prints out messages to confirm the operations applied.
            
        Returns
        -------
        Nothing
        """
        if self.status['ttx']:
            gnabar = 0.0
        else:
            gnabar = nstomho(gbar, self.somaarea)
        nach = self.status['na']
        if nach == 'nacncoop':
            soma().nacncoop.gbar = gnabar
            soma().nacncoop.KJ = 2000.
            soma().nacncoop.p = 0.25
            soma.ena = self.e_na
            if debug:
                print 'nacncoop gbar: ', soma().nacncoop.gbar
        elif nach == 'jsrna':
            soma().jsrna.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'jsrna gbar: ', soma().jsrna.gbar
        elif nach == 'nav11':
            soma().nav11.gbar = gnabar * 0.5
            soma.ena = self.e_na
            soma().nav11.vsna = 4.3
            if debug:
                print "Tuberculoventral using inva11"
            print 'nav11 gbar: ', soma().nav11.gbar
        elif nach == 'na':
            soma().na.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'na gbar: ', soma().na.gbar
        elif  nach == 'nacn':
            soma().nacn.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'nacn gbar: ', soma().nacn.gbar
        else:
            raise ValueError("Tuberculoventral setting Na channels: channel %s not known" % nach)

    # def add_axon(self):
    #     Cell.add_axon(self, self.soma, self.somaarea, self.c_m, self.R_a, self.axonsf)
    #
    # def add_dendrites(self):
    #     """
    #     Add 2 simple unbranched dendrites the basic model.
    #     The dendrites have some kht and ih current
    #     """
    #     print 'adding dendrites'
    #     cs = False  # not implemented outside here - internal Cesium.
    #     nDend = range(2) # these will be simple, unbranced, N=2 dendrites
    #     dendrites=[]
    #     for i in nDend:
    #         dendrites.append(h.Section(cell=self.soma))
    #     for i in nDend:
    #         dendrites[i].connect(self.soma)
    #         dendrites[i].L = 200. # length of the dendrite (not tapered)
    #         dendrites[i].diam = 1.0 # dendrite diameter
    #         dendrites[i].nseg = 21 # # segments in dendrites
    #         dendrites[i].Ra = 150 # ohm.cm
    #         dendrites[i].insert('kht')
    #         if cs is False:
    #             dendrites[i]().kht.gbar = 0.000 # a little Ht
    #         else:
    #             dendrites[i]().kht.gbar = 0.0
    #         dendrites[i].insert('leak') # leak
    #         dendrites[i]().leak.gbar = 0.001
    #         dendrites[i]().leak.erev = -78.0
    #         dendrites[i].insert('ihvcn') # some H current
    #         dendrites[i]().ihvcn.gbar = 0.0000
    #         dendrites[i]().ihvcn.eh = -43.0
    #         dendrites[i].ek = self.e_k
    #     self.maindend = dendrites
    #     self.status['dendrites'] = True
    #     self.add_section(self.maindend, 'maindend')



