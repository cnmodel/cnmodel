from __future__ import print_function
from neuron import h

from .cell import Cell
# from .. import synapses
from ..util import nstomho
from ..util import Params
import numpy as np
from .. import data

__all__ = ['MSO']


class MSO(Cell):
    
    celltype = 'mso'
    scaled = False
    
    @classmethod
    def create(cls, model='MSO-principal', **kwds):
        if model == 'MSO-principal':
            return MSOPrincipal(**kwds)
        else:
            raise ValueError ('MSO cell model %s is unknown', model)

    def make_psd(self, terminal, psd_type, **kwds):
        """
        Connect a presynaptic terminal to one post section at the specified location, with the fraction
        of the "standard" conductance determined by gbar.
        The default condition is designed to pass the unit test (loc=0.5)
        
        Parameters
        ----------
        terminal : Presynaptic terminal (NEURON object)
        
        psd_type : either simple or multisite PSD for MSO cell
        
        kwds: dictionary of options. 
            Two are currently handled:
            postsite : expect a list consisting of [sectionno, location (float)]
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
            if terminal.cell.celltype == 'bushy':
                # Max conductances for the glu mechanisms are calibrated by 
                # running `synapses/tests/test_psd.py`. The test should fail
                # if these values are incorrect
                self.AMPAR_gmax = data.get('bushy_synapse', species=self.species,
                        post_type=self.celltype, field='AMPAR_gmax')*1e3
                self.NMDAR_gmax = data.get('bushy_synapse', species=self.species,
                        post_type=self.celltype, field='NMDAR_gmax')*1e3
                self.Pr = data.get('bushy_synapse', species=self.species,
                        post_type=self.celltype, field='Pr')
                # adjust gmax to correct for initial Pr
                self.AMPAR_gmax = self.AMPAR_gmax/self.Pr
                self.NMDAR_gmax = self.NMDAR_gmax/self.Pr
                if 'AMPAScale' in kwds:  # normally, this should not be done!
                    self.AMPAR_gmax = self.AMPAR_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    self.NMDAR_gmax = self.NMDAR_gmax * kwds['NMDAScale']  # and NMDA... 
                return self.make_glu_psd(post_sec, terminal, self.AMPAR_gmax, self.NMDAR_gmax, loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.celltype, self.celltype))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)


class MSOPrincipal(MSO):
    """
    VCN MSO cell models.
        Using Rothman and Manis, 2003abc (Type II)
        MSO principal cell type
    """

    def __init__(self, morphology=None, decorator=None, nach=None,
                 ttx=False, species='guineapig', modelType=None, modelName=None, debug=False, temperature=None):
        """
        Create a MSO principal cell, using the default parameters for guinea pig from
        R&M2003, as a type II cell.
        Additional modifications to the cell can be made by calling methods below.
        
        Parameters
        ----------
        morphology : string (default: None)
            Name of a .hoc file representing the morphology. This file is used to constructe
            an electrotonic (cable) model. 
            If None (default), then a "point" (really, single cylinder) model is made, exactly according to RM03.
            
        decorator : Python function (default: None)
            decorator is a function that "decorates" the morphology with ion channels according
            to a set of rules.
            If None, a default set of channels is inserted into the first soma section, and the
            rest of the structure is "bare".
        
        nach : string (default: None)
            nach selects the type of sodium channel that will be used in the model. A channel mechanism
            by that name must exist. The default channel is set to 'nacn' (R&M03)
        
        temperature : float (default: 22)
            temperature to run the cell at. 
                 
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            This flag duplicates the effects of tetrodotoxin in the model. Currently, the flag is not implemented.
        
        temperature : float (default: None, sets to model default of 22)
            temperature (deg C) to run the cell at. Must be a valid temperature for the model.
        
        species: string (default 'guineapig')
            species defines the pattern of ion channel densities that will be inserted, according to 
            prior measurements in various species. Note that
            if a decorator function is specified, this argument is ignored as the decorator will
            specify the channel density.
            
        modelType: string (default: None)
            modelType specifies the subtype of the cell model that will be used (e.g., "II", "II-I", etc).
            modelType is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
            
        debug: boolean (default: False)
            When True, there will be multiple printouts of progress and parameters.
            
        Returns
        -------
            Nothing
                 
        """
        super(MSO, self).__init__()
        if species == 'guineapig':
            modelType = 'MSO-principal'
            modelName = 'MSO'
            dataset = 'MSO_principal_channels'
        else:
            raise ValueError(f"Species {species:s} not recognized for {self.celltype:s} cells")
        
        self.status = {self.somaname: True, 'axon': False, 'dendrites': False, 'pumps': False, 'hillock': False, 
                       'initialsegment': False, 'myelinatedaxon': False, 'unmyelinatedaxon': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'MSO',
                       'morphology': morphology, 'decorator': decorator, 'temperature': temperature}

        self.debug = debug
        self.spike_threshold = -40
        soma = self.do_morphology(morphology)
        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.natype
        
        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            self.mechanisms = ['klt', 'kht', 'ihvcn', 'leak', self.pars.natype]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ena = self.e_na
            self.soma.ek = self.e_k
            self.soma().ihvcn.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.c_m = 0.9
            self.species_scaling(silent=True)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)

        if debug:
            print("   << Created cell >>")

    def get_cellpars(self, dataset, species='guineapig', modelType='principal'):
        cellcap = data.get(dataset, species=species, model_type=modelType,
            field='soma_Cap')
        chtype = data.get(dataset, species=species, model_type=modelType,
            field='soma_na_type')
        pars = Params(cap=cellcap, natype=chtype)
        for g in ['soma_kht_gbar', 'soma_klt_gbar', 'soma_ih_gbar', 'soma_leak_gbar']:
            pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
            field=g))
        return pars

    def species_scaling(self, silent=True):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        This scaling should be used ONLY for point models, as no other compartments
        are scaled.
        
        This scaling routine also sets the temperature for the model to a default value. Some models
        can be run at multiple temperatures, and so a default from one of the temperatures is used.
        The calling cell.set_temperature(newtemp) will change the conductances and reinitialize
        the cell to the new temperature settings.
        
        Parameters
        ----------
        species : string (default: 'guineapig')
            name of the species to use for scaling the conductances in the base point model
            Must be one of mouse, cat, guineapig
        
        modelType: string (default: 'principal')
            definition of model type from RM03 models, principal cell for mso 
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        
        """
        assert self.scaled is False  # block double scaling!
        self.scaled = True
            
        soma = self.soma
            
        if self.status['species'] == 'guineapig':
            print(f"Setting conductances for guinea pig {self.status['modelType']:s} MSO cell, based on Rothman and Manis, 2003 bushy cell")
            self._valid_temperatures = (22., 38.)
            if self.status['temperature'] is None:
                self.status['temperature'] = 22. 
            self.i_test_range = {'pulse': (-0.4, 0.4, 0.02)}
            self.vrange = [-70., -55.]  # set a default vrange for searching for rmp
            sf = 1.0
            if self.status['temperature'] == 38.:  # adjust for 2003 model conductance levels at 38
                sf = 2  # Q10 of 2, 22->38C. (p3106, R&M2003c)
                # note that kinetics are scaled in the mod file.
            self.set_soma_size_from_Cm(self.pars.cap)

            self.adjust_na_chans(soma, sf=sf)
            soma().kht.gbar = nstomho(self.pars.soma_kht_gbar, self.somaarea)
            soma().klt.gbar = nstomho(self.pars.soma_klt_gbar, self.somaarea)
            soma().ihvcn.gbar = nstomho(self.pars.soma_ih_gbar, self.somaarea)
            soma().leak.gbar = nstomho(self.pars.soma_leak_gbar, self.somaarea)

            self.axonsf = 0.57
            
        else:
            raise ValueError(f"Species {self.status['species']:s} or species-type {self.status['modelType']:s} is not recognized for MSO cells")
        self.check_temperature()
#        self.cell_initialize(vrange=self.vrange)  # no need to do this just yet.
        if not silent:
           print(' set cell as: ', {self.status['species']:s})

