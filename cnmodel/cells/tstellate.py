from __future__ import print_function
from neuron import h
import numpy as np

from .cell import Cell
#from .. import synapses
from ..util import nstomho
from ..util import Params
from .. import data

__all__ = ['TStellate', 'TStellateRothman'] 

"""
Class to instantiate T stellate (planar multipolar) cells from the cochlear nucleus.

02 May 2019: Old Tstellate Nav11 removed (can instantiate by selecting nach during creation instead; data
is read from table, so nav11 channel is used for mouse I-t model with XM13_channels)


"""

class TStellate(Cell):
    
    celltype = 'tstellate'
    scaled = False
    
    @classmethod
    def create(cls, model='RM03', **kwds):
        if model == 'RM03':  # original Rothman-Manis 2003, 22C, point cell, extendable
            return TStellateRothman(**kwds)
        else:
            raise ValueError ('TStellate type %s is unknown', type)
        

    def make_psd(self, terminal, psd_type, **kwds):
        """
        Connect a presynaptic terminal to one post section at the specified location, with the fraction
        of the "standard" conductance determined by gbar.
        The default condition is to try to pass the default unit test (loc=0.5)
        
        Scaling is corrected by initial release probability now.
        
        Parameters
        ----------
        terminal : Presynaptic terminal (NEURON object)
        
        psd_type : either simple or multisite PSD for bushy cell
        
        kwds: dict of options.
            Available options:
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
        #print('cells/tstellaty.py psd type: ', psd_type)
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
                # old values:
                # AMPA_gmax = 0.22479596944138733*1e3  # factor of 1e3 scales to pS (.mod mechanisms) from nS.
                # NMDA_gmax = 0.12281291946623739*1e3
                if 'AMPAScale' in kwds:
                    self.AMPAR_gmax = self.AMPAR_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    self.NMDAR_gmax = self.NMDAR_gmax*kwds['NMDAScale']
                return self.make_glu_psd(post_sec, terminal, self.AMPAR_gmax, self.NMDAR_gmax, loc=loc)
            elif terminal.cell.celltype == 'dstellate':
                self.ds_gmax = data.get('dstellate_synapse', species=self.species,
                        post_type=self.celltype, field='gly_gmax')*1e3
                return self.make_gly_psd(post_sec, terminal, psdtype='glyfast', loc=loc, gmax=self.ds_gmax)
            elif terminal.cell.celltype == 'tuberculoventral':
                self.tv_gmax = data.get('tuberculoventral_synapse', species=self.species,
                        post_type=self.celltype, field='gly_gmax')*1e3
                return self.make_gly_psd(post_sec, terminal, psdtype='glyfast', loc=loc, gmax=self.tv_gmax)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.celltype, self.celltype))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)


class TStellateRothman(TStellate):
    """
    VCN T-stellate base model.
    Rothman and Manis, 2003abc (Type I-c, Type I-t)
    """
    def __init__(self, morphology=None, decorator=None, nach=None,
                ttx=False, temperature=None,
                species='guineapig', modelType=None, modelName=None, 
                debug=False):
        """
        Initialize a planar stellate (T-stellate) cell, using the default parameters for guinea pig from
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
    
        nach : string (default: None)
            nach selects the type of sodium channel that will be used in the model. A channel mechanism
            by that name must exist. The default is 'nacn', from R&M2003.
    
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
    
        species: string (default 'guineapig')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored.

        modelName: string (default: None)
            modelName specifies the source conductance pattern (RM03, XM13, etc).
            modelName is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
                                     
        modelType: string (default: None)
            modelType specifies the type of the model that will be used (e.g., "I-c", "I-t").
            modelType is passed to the decorator, or to species_scaling to adjust point models.
        
        debug: boolean (default: False)
            debug is a boolean flag. When set, there will be multiple printouts of progress and parameters.
        
        Returns
        -------
        Nothing
        """
        
        super(TStellateRothman, self).__init__()
        if modelType is None:
            modelType = 'I-c'
        
        if species == 'guineapig':
            modelName = 'RM03'
            temp = 22.
            dataset = 'RM03_channels'
        
        elif species == 'mouse':
            temp = 34.
            if modelName is None:
                modelName = 'XM13'
            if modelName == 'XM13':
                dataset = 'XM13_channels'
            elif modelName  == 'XM13nacncoop':
                dataset = 'XM13_channels_nacncoop'
            else:
                raise ValueError(f"ModelName {self.status['modelName']:s} not recognized for {self.celltype:s} cells")
        else:
            raise ValueError(f"Species {species:s} not recognized for {self.celltype:s} cells")

        self.debug = debug
        self.status = {'species': species, 'cellClass': self.celltype, 'modelType': modelType, 'modelName': modelName,
                       self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'ttx': ttx, 'name': self.celltype,
                       'morphology': morphology, 'decorator': decorator, 'temperature': None}
        self.c_m = 0.9  # default in units of uF/cm^2
        self.spike_threshold = -40.  # matches threshold in released CNModel (set in base cell class)

        if self.debug:
            print( 'model type, model name, species: ', modelType, modelName, species, nach)

        self._valid_temperatures = (temp, )
        if self.status['temperature'] == None:
            self.status['temperature'] = temp

        soma = self.do_morphology(morphology)

        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.natype
        
        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['kht', 'ka', 'ihvcn', 'leak', self.pars.natype]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ek = self.e_k
            self.soma.ena = self.e_na
            self.soma().ihvcn.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.species_scaling(silent=True)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()

        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)
        if self.debug:
                print("<< T-stellate: JSR Stellate Type 1 cell model created >>")

    def get_cellpars(self, dataset, species='guineapig', modelType='I-c'):
        """
        Retrieve parameters for the specifed model type and species from the data tables
        
        dataset : str (no default)
            name of the data table to use
        
        species : str (default: 'guineapig')
            Species table to use
        
        modelType : str (default: 'I-c')
            Model type to get parameters from the table.
        
        """
        cellcap = data.get(dataset, species=species, model_type=modelType,
            field='soma_Cap')
        chtype = data.get(dataset, species=species, model_type=modelType,
            field='na_type')
        pars = Params(cap=cellcap, natype=chtype)

        if self.status['modelName'] == 'RM03':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'ka_gbar', 'ih_gbar', 'leak_gbar', 'leak_erev', 'ih_eh', 'e_k', 'e_na']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        elif self.status['modelName'] == 'XM13':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'ka_gbar', 'ihvcn_gbar', 'leak_gbar', 'leak_erev', 'ih_eh', 'e_k', 'e_na']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        # elif self.status['modelName'] == 'mGBC':
        #     for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'ka_gbar', 'ihvcn_gbar', 'leak_gbar', 'leak_erev', 'ih_eh', 'e_k', 'e_na']:
        #         pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
        #             field=g))
        else:
            raise ValueError(f"get_cellpars: Model name {self.status['modelName']} is not yet implemented for cell type {self.celltype.title():s}")
            
        if self.debug:
            pars.show()
        return pars
    

    def species_scaling(self, silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
        
        This scaling routine also sets the temperature for the model to a default value. Some models
        can be run at multiple temperatures, and so a default from one of the temperatures is used.
        The calling cell.set_temperature(newtemp) will change the conductances and reinitialize
        the cell to the new temperature settings.
        
        Parameters
        ----------
        species : string (default: 'guineapig')
            name of the species to use for scaling the conductances in the base point model
            Must be one of mouse, cat, guineapig
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        """
        
        assert self.scaled is False  # block double scaling!
        self.scaled = True
        
        soma = self.soma

        if self.status['species'] == 'mouse': #  and modelType == 'I-c':
            # use conductance levels from Cao et al.,  J. Neurophys., 2007.
            # model description in Xie and Manis 2013. Note that
            # conductances were not scaled for temperature (rates were)
            # so here we reset the default Q10's for conductance (g) to 1.0
            if self.status['modelType'] not in ['I-c', 'I-t']:
                raise ValueError(f"\nModel type {self.status['modelType']:s} is not implemented for mouse {self.celltype.title():s}")
            if self.debug:
                print(f"  Setting Conductances for mouse {self.status['modelType']:s} {self.celltype.title():s} cell, (modified from Xie and Manis, 2013)")
            # self.c_m = 0.9  # default in units of F/cm^2
            self.vrange = [-75., -55.]
            self.i_test_range={'pulse': (-1.0, 1.0, 0.05)}
            self._valid_temperatures = (34.,)
            if self.status['temperature'] is None:
                self.set_temperature(34.)

            self.set_soma_size_from_Cm(self.pars.cap)
            self.adjust_na_chans(soma, sf=1.0)
            soma().kht.gbar = nstomho(self.pars.kht_gbar, self.somaarea)
            soma().ka.gbar = nstomho(self.pars.ka_gbar, self.somaarea)
            soma().ihvcn.gbar = nstomho(self.pars.ihvcn_gbar, self.somaarea)
            soma().ihvcn.eh = self.pars.ih_eh # Rodrigues and Oertel, 2006
            soma().leak.gbar = nstomho(self.pars.leak_gbar, self.somaarea)
            soma().leak.erev = self.pars.leak_erev
            self.e_k = self.pars.e_k
            self.e_na = self.pars.e_na
            soma.ena = self.pars.e_na
            soma.ek = self.e_k
            self.axonsf = 0.5
            
        elif self.status['species'] == 'guineapig':
            if self.status['modelType'] not in ['I-c', 'I-t']:
                raise ValueError(f"\nModel type {self.status['modelType']:s} is not implemented for mouse {self.celltype.title():s} cells")
            if self.debug:
                print(f"  Setting Conductances for Guinea Pig {self.status['modelType']:s} {self.celltype.title():s}  cell, Rothman and Manis, 2003")
            self.c_m = 0.9  # default in units of F/cm^2
            self.vrange = [-75., -55.]
            self.i_test_range={'pulse': (-0.15, 0.15, 0.01)}
            self._valid_temperatures = (22., 38.)
            if self.status['temperature'] is None:
                self.set_temperature(22.)
            sf = 1.0
            if self.status['temperature'] == 38.:  # adjust for 2003 model conductance levels at 38
                sf = 3.03  # Q10 of 2, 22->38C. (p3106, R&M2003c)
                # note that kinetics are scaled in the mod file.
            self.set_soma_size_from_Cm(self.pars.cap)
            self.adjust_na_chans(soma, sf=sf)
            soma().kht.gbar = nstomho(self.pars.kht_gbar, self.somaarea)
            soma().ka.gbar = nstomho(self.pars.ka_gbar, self.somaarea)
            soma().ihvcn.gbar = nstomho(self.pars.ih_gbar, self.somaarea)
            soma().leak.gbar = nstomho(self.pars.leak_gbar, self.somaarea)
            soma().leak.erev = self.pars.leak_erev
            self.axonsf = 0.5

        else:
            raise ValueError(f"Species {self.status['species']:s} or species-type {self.status['modelType']:s} is not recognized for T-stellate cells")

        self.check_temperature()

    def get_distancemap(self):
        """
        This structure should be in a data table
        """
        return {'dend': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 200.},
                         'kht': {'gradient': 'llinear', 'gminf': 0., 'lambda': 200.},
                         'nav11': {'gradient': 'linear', 'gminf': 0., 'lambda': 200.}}, # linear with distance, gminf (factor) is multiplied by gbar
                'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                         'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                         'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 200.}}, # gradients are: flat, linear, exponential
                }

    def add_axon(self):
        Cell.add_axon(self, self.soma, self.somaarea, self.c_m, self.R_a, self.axonsf)

    def add_dendrites(self):
        """
        Add simple unbranched dendrites to basic Rothman Type I models.
        The dendrites have some kht and ih current
        """
        cs = False  # not implemented outside here - internal Cesium.
        nDend = range(4) # these will be simple, unbranced, N=4 dendrites
        dendrites=[]
        for i in nDend:
            dendrites.append(h.Section(cell=self.soma))
        for i in nDend:
            dendrites[i].connect(self.soma)
            dendrites[i].L = 200 # length of the dendrite (not tapered)
            dendrites[i].diam = 1.5 # dendrite diameter
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

