from __future__ import print_function
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
    
    celltype = 'tuberculoventral'
    scaled = False
    
    @classmethod
    def create(cls, model='TVmouse', **kwds):
        if model in ['TVmouse', 'I']:
            return Tuberculoventral(**kwds)
        elif model == 'dummy':
            return DummyTuberculoventral(**kwds)
        else:
            raise ValueError ('Tuberculoventral type %s is unknown', model)

    def __init__(self):
        Cell.__init__(self)
        self.spike_source = None  # used by DummyTuberculoventral to connect VecStim to terminal

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
            if terminal.cell.celltype in ['sgc', 'dstellate', 'tuberculoventral']:
                weight = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='weight')
                tau1 = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='tau1')
                tau2 = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='tau2')
                erev = data.get('%s_synapse' % terminal.cell.celltype, species=self.species,
                        post_type=self.celltype, field='erev')
                return self.make_exp2_psd(post_sec, terminal, weight=weight, loc=loc,
                        tau1=tau1, tau2=tau2, erev=erev)
            else:
                raise TypeError("Cannot make simple PSD for %s => %s" % 
                            (terminal.cell.celltype, self.celltype))

        elif psd_type == 'multisite':
            if terminal.cell.celltype == 'sgc':
                # Max conductances for the glu mechanisms are calibrated by 
                # running `synapses/tests/test_psd.py`. The test should fail
                # if these values are incorrect
                self.AMPAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.celltype, field='AMPAR_gmax')*1e3
                self.NMDAR_gmax = data.get('sgc_synapse', species=self.species,
                        post_type=self.celltype, field='NMDAR_gmax')*1e3
                self.Pr = data.get('sgc_synapse', species=self.species,
                        post_type=self.celltype, field='Pr')
                # adjust gmax to correct for initial Pr
                self.AMPAR_gmax = self.AMPAR_gmax/self.Pr
                self.NMDAR_gmax = self.NMDAR_gmax/self.Pr
                if 'AMPAScale' in kwds:
                    self.AMPA_gmax = self.AMPA_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    self.NMDA_gmax = self.NMDA_gmax*kwds['NMDAScale']
                return self.make_glu_psd(post_sec, terminal, self.AMPAR_gmax, self.NMDAR_gmax, loc=loc)
            elif terminal.cell.celltype == 'dstellate':  # WBI input -Voigt, Nelken, Young
                return self.make_gly_psd(post_sec, terminal, psdtype='glyfast', loc=loc)
            elif terminal.cell.celltype == 'tuberculoventral':  # TV cells talk to each other-Kuo et al.
                return self.make_gly_psd(post_sec, terminal, psdtype='glyfast', loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.celltype, self.celltype))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)

    def make_terminal(self, post_cell, term_type, **kwds):
        pre_sec = self.soma
        if term_type == 'simple':
            return synapses.SimpleTerminal(pre_sec, post_cell, spike_source=self.spike_source, 
                                            **kwds)
        elif term_type == 'multisite':
            if post_cell.celltype in ['dstellate', 'tuberculoventral', 'pyramidal', 'bushy', 'tstellate']:
                nzones = data.get('tuberculoventral_synapse', species=self.species,
                        post_type=post_cell.celltype, field='n_rsites')
                delay = data.get('tuberculoventral_synapse', species=self.species,
                        post_type=post_cell.celltype, field='delay')
            else:
                raise NotImplementedError("No knowledge as to how to connect tuberculoventral cell to cell type %s" %
                                        type(post_cell))
            pre_sec = self.soma
            return synapses.StochasticTerminal(pre_sec, post_cell, nzones=nzones, spike_source=self.spike_source,
                                            delay=delay, **kwds)
        else:
            raise ValueError("Unsupported terminal type %s" % term_type)


class Tuberculoventral(Tuberculoventral):
    """
    Tuberculoventral Neuron (DCN) base model
    Adapted from T-stellate model, using target parameters from Kuo et al. J. Neurophys. 2012
    """
    def __init__(self, morphology=None, decorator=None, nach=None, ttx=False,
                species='mouse', modelType=None, modelName=None, debug=False):
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
        if species == 'mouse':
            temp = 34.
            if modelName is None:
                modelName = 'TVmouse'
            if modelName == 'TVmouse':
                dataset = 'TV_channels'
            else:
                raise ValueError(f"ModelName {self.status['modelName']:s} not recognized for mouse T-stellate cells")
        else:
            raise ValueError(f"Species {species:s} not recognized for {self.celltype:s} cells")
        self.debug = debug  
        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'modelName': modelName,
                       'ttx': ttx, 'name': 'Tuberculoventral',
                       'morphology': morphology, 'decorator': decorator, 'temperature': None}

        self.do_morphology(morphology)
        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.soma_natype

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['kht', 'ka', 'ihvcn', 'leak', self.pars.soma_natype]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.species_scaling(silent=True)  # adjust the default parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)
        if self.debug:
                print("<< Mouse TV cell created>>")

    def get_cellpars(self, dataset, species='mouse', modelType='TVmouse'):
        cellcap = data.get(dataset, species=species, model_type=modelType,
            field='soma_Cap')
        chtype = data.get(dataset, species=species, model_type=modelType,
            field='soma_na_type')
        pars = Params(soma_cap=cellcap, soma_natype=chtype)
        for g in ['soma_nacncoop_gbar', 'soma_kht_gbar', 'soma_ka_gbar',
                  'soma_ihvcn_gbar', 'soma_ihvcn_eh',
                  'soma_leak_gbar', 'soma_leak_erev',
                  'soma_e_k', 'soma_e_na']:
            pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
            field=g))
        return pars
        
    def species_scaling(self, silent=True):
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
        assert self.scaled is False  # block double scaling!
        self.scaled = True
        
        soma = self.soma

        if self.status['species'] == 'mouse':
            if self.status['modelType'] not in ['TVmouse']:
                raise ValueError('\nModel type %s is not implemented for mouse Tstellate cells' % self.status['modelType'])
            if self.debug:
                print(f"  Setting Conductances for mouse {self.status['modelType']:s} Tstellate cell, (modified from Xie and Manis, 2013)")
                
            """#From Kuo 150 Mohm, 10 msec tau
            Firing at 600 pA about 400 Hz
            These values from brute_force runs, getting 380 Hz at 600 pA at 35C
            Input resistance and vm is ok, time constnat is short
                *** Rin:       168  tau:       7.8   v:  -68.4
            Attempts to get longer time constant - cannot keep rate up.
            """
            # Adapted from TStellate model type I-c'
            self.i_test_range = {'pulse': [(-0.35, 1.0, 0.05), (-0.04, 0.01, 0.01)]}
            self.vrange = [-80., -58.]  # set a default vrange for searching for rmp
            self._valid_temperatures = (34.,)
            if self.status['temperature'] is None:
                self.set_temperature(34.)

            self.set_soma_size_from_Cm(self.pars.soma_cap)
            self.adjust_na_chans(soma)
            soma().kht.gbar = nstomho(self.pars.soma_kht_gbar, self.somaarea)
            soma().ka.gbar = nstomho(self.pars.soma_ka_gbar, self.somaarea)
            soma().ihvcn.gbar = nstomho(self.pars.soma_ihvcn_gbar, self.somaarea)
            soma().ihvcn.eh = self.pars.soma_ihvcn_eh
            soma().leak.gbar = nstomho(self.pars.soma_leak_gbar, self.somaarea)
            soma().leak.erev = self.pars.soma_leak_erev
            self.e_leak = self.pars.soma_leak_erev
            self.soma.ek = self.e_k = self.pars.soma_e_k
            self.soma.ena = self.e_na = self.pars.soma_e_na

            self.axonsf = 0.5
        else:
            raise ValueError(f"Species {self.status['species']:s} or species modeltype {self.status['modelType']:s} is not recognized for {self.celltype:s} cells")


        self.check_temperature()


#     def channel_manager(self, modelType='TVmouse'):
#         """
#         This routine defines channel density maps and distance map patterns
#         for each type of compartment in the cell. The maps
#         are used by the ChannelDecorator class (specifically, it's private
#         _biophys function) to decorate the cell membrane.
#
#         Parameters
#         ----------
#         modelType : string (default: 'RM03')
#             A string that defines the type of the model. Currently, 3 types are implemented:
#             RM03: Rothman and Manis, 2003 somatic densities for guinea pig
#             XM13: Xie and Manis, 2013, somatic densities for mouse
#             XM13PasDend: XM13, but with only passive dendrites, no channels.
#
#         Returns
#         -------
#         Nothing
#
#         Notes
#         -----
#
#         This routine defines the following variables for the class:
#
#             - conductances (gBar)
#             - a channelMap (dictonary of channel densities in defined anatomical compartments)
#             - a current injection range for IV's (when testing)
#             - a distance map, which defines how selected conductances in selected compartments
#                 will change with distance. This includes both linear and exponential gradients,
#                 the minimum conductance at the end of the gradient, and the space constant or
#                 slope for the gradient.
#
#         """
#         if modelType == 'TVmouse':
#             print('decorate as tvmouse')
# #            totcap = 95.0E-12  # Tuberculoventral cell (type I), based on stellate, adjusted for Kuo et al. TV firing
#             self.set_soma_size_from_Section(self.soma)
#             totcap = self.totcap
#             refarea = self.somaarea # totcap / self.c_m  # see above for units
#             self.gBar = Params(nabar=1520.0E-9/refarea,
#                                khtbar=160.0E-9/refarea,
#                                kltbar=0.0E-9/refarea,
#                                kabar=65.0/refarea,
#                                ihbar=1.25E-9/refarea,
#                                leakbar=5.5E-9/refarea,
#             )
#             self.channelMap = {
#                 'axon': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar,
#                          'ihvcn': 0., 'leak': self.gBar.leakbar / 4.},
#                 'hillock': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
#                              'ihvcn': 0., 'leak': self.gBar.leakbar, },
#                 'initseg': {'nacn': self.gBar.nabar, 'klt': 0., 'kht': self.gBar.khtbar,
#                             'ihvcn': self.gBar.ihbar / 2.,
#                             'leak': self.gBar.leakbar, },
#                 self.somaname: {'nacn': self.gBar.nabar, 'klt': self.gBar.kltbar,
#                          'kht': self.gBar.khtbar, 'ihvcn': self.gBar.ihbar,
#                          'leak': self.gBar.leakbar, },
#                 'dend': {'nacn': self.gBar.nabar / 2.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.5,
#                          'ihvcn': self.gBar.ihbar / 3., 'leak': self.gBar.leakbar * 0.5, },
#                 'apic': {'nacn': 0.0, 'klt': 0., 'kht': self.gBar.khtbar * 0.2,
#                          'ihvcn': self.gBar.ihbar / 4.,
#                          'leak': self.gBar.leakbar * 0.2, },
#             }
#             self.irange = np.linspace(-0.3, 0.6, 10)
#             self.distMap = {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
#                                      'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
#                             'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
#                                      'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # gradients are: flat, linear, exponential
#                             }
#         else:
#             raise ValueError('model type %s is not implemented' % modelType)


class DummyTuberculoventral(Tuberculoventral):
    """ Tuberculoventral cell class with no cell body; this cell only replays a predetermined
    spike train. Useful for testing, or replacing spike trains to determine
    the importance of spike structures within a network.
    """
    def __init__(self, cf=None, species='mouse'):
        """
        Parameters
        ----------
        cf : float (default: None)
            Required: the characteristic frequency for the TV cell
            Really just for reference.

        """

        Tuberculoventral.__init__(self)
        self.vecstim = h.VecStim()
        
        # this causes the terminal to receive events from the VecStim:
        self.spike_source = self.vecstim
        
        # just an empty section for holding the terminal
        self.add_section(h.Section(), self.somaname)
        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': None, 'species': species, 'modelType': 'Dummy', 'modelName': 'DummyTuberculoventral',
                       'ttx': None, 'name': 'DummyTuberculoventral',
                       'morphology': None, 'decorator': None, 'temperature': None}
        print("<< Tuberculoventral: Dummy Tuberculoventral Cell created >>")
        

    def set_spiketrain(self, times):
        """ Set the times of spikes (in seconds) to be replayed by the cell.
        """
        self._spiketrain = times
        self._stvec = h.Vector(times)
        self.vecstim.play(self._stvec)


