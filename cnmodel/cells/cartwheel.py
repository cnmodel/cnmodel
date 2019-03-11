from __future__ import print_function
from neuron import h
from ..util import nstomho
from ..util import Params
import numpy as np
from .cell import Cell
from .. import data
from .. import synapses

__all__ = ['Cartwheel', 'CartwheelDefault']

class Cartwheel(Cell):

    @classmethod
    def create(cls, model='CW', **kwds):
        if model == 'CW':
            return CartwheelDefault(**kwds)
        else:
            raise ValueError ('Carthweel model is unknown', model)

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
        self.pre_sec = terminal.section
        pre_cell = terminal.cell
        post_sec = self.soma

        
        if psd_type == 'simple':
            if terminal.cell.type in ['sgc', 'dstellate', 'tuberculoventral', 'cartwheel']:
                weight = data.get('%s_synapse' % terminal.cell.type, species=self.species,
                        post_type=self.type, field='weight')
                tau1 = data.get('%s_synapse' % terminal.cell.type, species=self.species,
                        post_type=self.type, field='tau1')
                tau2 = data.get('%s_synapse' % terminal.cell.type, species=self.species,
                        post_type=self.type, field='tau2')
                erev = data.get('%s_synapse' % terminal.cell.type, species=self.species,
                        post_type=self.type, field='erev')
                return self.make_exp2_psd(post_sec, terminal, weight=weight, loc=loc,
                        tau1=tau1, tau2=tau2, erev=erev)
            else:
                raise TypeError("Cannot make simple PSD for %s => %s" % 
                            (terminal.cell.type, self.type))

        else:
            raise ValueError("Unsupported psd type %s for cartwheel cell (inputs not implemented yet)" % psd_type)

    def make_terminal(self, post_cell, term_type, **kwds):
        if term_type == 'simple':
            return synapses.SimpleTerminal(self.pre_sec, post_cell, 
                                            **kwds)
        elif term_type == 'multisite':
            if post_cell.type in ['tuberculoventral', 'pyramidal']:
                nzones = data.get('cartwheel_synapse', species=self.species,
                        post_type=post_cell.type, field='n_rsites')
                delay = data.get('cartwheel_synapse', species=self.species,
                        post_type=post_cell.type, field='delay')
            else:
                raise NotImplementedError("No knowledge as to how to connect cartwheel cell to cell type %s" %
                                        type(post_cell))
            pre_sec = self.soma
            return synapses.StochasticTerminal(pre_sec, post_cell, nzones=nzones, spike_source=self.spike_source,
                                            delay=delay, **kwds)
        else:
            raise ValueError("Unsupported terminal type %s" % term_type)

class CartwheelDefault(Cartwheel, Cell):
    """
    DCN cartwheel cell model.
    
    """
    def __init__(self, morphology=None, decorator=None, ttx=False, nach=None,
                 species='mouse', modelType=None, debug=False):
        """        
        Create cartwheel cell model, based on a Purkinje cell model from Raman.
        There are no variations available for this model.
        
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
            by that name must exist. The default is naRsg, a resurgent sodium channel model.
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            This flag duplicates the effects of tetrodotoxin in the model. Currently, the flag is not implemented.
        
        species: string (default 'rat')
            species defines the pattern of ion channel densities that will be inserted, according to 
            prior measurements in various species. Note that
            if a decorator function is specified, this argument is ignored as the decorator will
            specify the channel density.
            
        modelType: string (default: None)
            modelType specifies the subtype of the cell model that will be used.
            modelType is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
            Only type "I" is recognized for the cartwheel cell model.
            
        debug: boolean (default: False)
            When True, there will be multiple printouts of progress and parameters.
            
        Returns
        -------
            Nothing
        """
        super(CartwheelDefault, self).__init__()
        if modelType == None:
            modelType = 'I'
        if nach == None:
            nach = 'naRsg'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'Cartwheel',
                       'morphology': morphology, 'decorator': decorator, 'temperature': None}

        self.i_test_range = {'pulse': (-0.2, 0.2, 0.02)}
       # self.spike_threshold = 0
        self.vrange = [-75., -52.]  # set a default vrange for searching for rmp

        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            soma = h.Section(name="Cartwheel_Soma_%x" % id(self)) # one compartment of about 29000 um2
#            cm = 1
            soma.nseg = 1
            self.add_section(soma, 'soma')
        else:
            """
            instantiate a structured model with the morphology as specified by 
            the morphology file
            """
            self.set_morphology(morphology_file=morphology)

        # decorate the morphology with ion channels
        if decorator is None:   # basic model, only on the soma
            # v_potassium = -80       # potassium reversal potential
            # v_sodium = 50           # sodium reversal potential

            self.mechanisms = ['naRsg', 'bkpkj', 'hpkj', 'kpkj', 'kpkj2',
                               'kpkjslow', 'kpksk', 'lkpkj', 'cap']
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.insert('cadiff')
            self.species_scaling(silent=True, species=species, modelType=modelType)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)
        
        if debug:
            print( "<< Cartwheel: Modified version of Raman Purkinje cell model created >>")

    def get_cellpars(self, dataset, species='guineapig', celltype='II'):
        somaDia = data.get(dataset, species=species, cell_type=celltype,
            field='soma_Dia')
        chtype = data.get(dataset, species=species, cell_type=celltype,
            field='soma_na_type')
        pcabar = data.get(dataset, species=species, cell_type=celltype,
            field='soma_pcabar')
        pars = Params(soma_Dia=somaDia, soma_natype=chtype, soma_pcabar=pcabar)
        for g in ['soma_narsg_gbar', 'soma_kpkj_gbar', 'soma_kpkj2_gbar', 'soma_kpkjslow_gbar',
                  'soma_kpksk_gbar', 'soma_lkpkj_gbar', 'soma_bkpkj_gbar', 'soma_hpkj_gbar',
                  'soma_hpkj_eh','soma_lkpkj_e', 'soma_e_k', 'soma_e_na', 'soma_e_ca',
                  ]:
            pars.additem(g,  data.get(dataset, species=species, cell_type=celltype,
            field=g))
        return pars
        
    def species_scaling(self, silent=True, species='mouse', modelType='I'):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        This scaling should be used ONLY for point models, as no other compartments
        are scaled.
        
        Parameters
        ----------
        species : string (default: 'rat')
            name of the species to use for scaling the conductances in the base point model
            Must be one of mouse, cat, guineapig
        
        modelType: string (default: 'I')
            definition of model type from RM03 models, type II or type II-I
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        
        Note
        ----
            For the cartwheel cell model, there is only a single scaling recognized. 
        """        
        if species is not 'mouse':
            raise ValueError ('Cartwheel species: only "mouse" is recognized')
        if modelType is not 'I':
            raise ValueError ('Cartwheel modelType: only "I" is recognized')
        self._valid_temperatures = (34.,)
        if self.status['temperature'] is None:
            self.set_temperature(34.)
        
        pars = self.get_cellpars('CW_channels', species=species, celltype='cartwheel')
        self.set_soma_size_from_Diam(pars.soma_Dia)
        self.soma().bkpkj.gbar = nstomho(pars.soma_bkpkj_gbar, self.somaarea)
        self.soma().hpkj.gbar = nstomho(pars.soma_hpkj_gbar, self.somaarea)
        self.soma().kpkj.gbar = nstomho(pars.soma_kpkj_gbar, self.somaarea)
        self.soma().kpkj2.gbar = nstomho(pars.soma_kpkj2_gbar, self.somaarea)
        self.soma().kpkjslow.gbar = nstomho(pars.soma_kpkjslow_gbar, self.somaarea)
        self.soma().kpksk.gbar = nstomho(pars.soma_kpksk_gbar, self.somaarea)
        self.soma().lkpkj.gbar = nstomho(pars.soma_lkpkj_gbar, self.somaarea)
        self.soma().naRsg.gbar = nstomho(pars.soma_narsg_gbar, self.somaarea)
        self.soma().cap.pcabar = pars.soma_pcabar
        self.soma().ena = pars.soma_e_na # 50
        self.soma().ek = pars.soma_e_k # -80
        self.soma().lkpkj.e = pars.soma_lkpkj_e  #-65
        self.soma().hpkj.eh = pars.soma_hpkj_eh  #-43
        self.soma().eca = pars.soma_e_ca # 50
        
        self.status['na'] = pars.soma_natype
        self.status['species'] = species
        self.status['modelType'] = modelType
        self.check_temperature()
        if not silent:
            print( 'set cell as: ', species)
            print( ' with Vm rest = %f' % self.vm0)
       # print 'set up'
        
    def i_currents(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point.
        Overrides i_currents in cells.py, because this model uses conductances
        that are not specified in the default cell mode.
        
        Parameters
        ----------
        V : float, mV (no default)
            Voltage at which the current for each conductance is computed.
        
        Returns
        -------
        I : float, nA
             The sum of the currents at steady-state for all of the conductances.
        """
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = V
        h.celsius = self.status['temperature']
        h.finitialize()
        self.ix = {}

        if 'naRsg' in self.mechanisms:
             self.ix['naRsg'] = self.soma().naRsg.gna*(V - self.soma().ena)
        if 'cap' in self.mechanisms:
            a = self.soma().cap.pcabar*self.soma().cap.minf
            self.ix['cap'] = a * self.ghk(V, self.soma().cao, self.soma().cai, 2)
        if 'kpkj' in self.mechanisms:
             self.ix['kpkj'] = self.soma().kpkj.gk*(V - self.soma().ek)
        if 'kpkj2' in self.mechanisms:
             self.ix['kpkj2'] = self.soma().kpkj2.gk*(V - self.soma().ek)
        if 'kpkjslow' in self.mechanisms:
             self.ix['kpkjslow'] = self.soma().kpkjslow.gk*(V - self.soma().ek)
        if 'kpksk' in self.mechanisms:
             self.ix['kpksk'] = self.soma().kpksk.gk*(V - self.soma().ek)
        if 'bkpkj' in self.mechanisms:
             self.ix['bkpkj'] = self.soma().bkpkj.gbkpkj*(V - self.soma().ek)
        if 'hpkj' in self.mechanisms:
             self.ix['hpkj'] = self.soma().hpkj.gh*(V - self.soma().hpkj.eh)
        # leak
        if 'lkpkj' in self.mechanisms:
            self.ix['lkpkj'] = self.soma().lkpkj.gbar*(V - self.soma().lkpkj.e)
        return np.sum([self.ix[i] for i in self.ix])

    def ghk(self, v, ci, co, z):
        """
        GHK flux equation, used to calculate current density through calcium channels
        rather than standard Nernst equation.
        
        Parameters
        ----------
        v : float, mV
            voltage for GHK calculation
        ci : float, mM
            internal ion concentration
        co : float, mM
            external ion concentraion
        z : float, no units
            valence
        
        Returns
        -------
        flux : A/m^2

        """
        F = 9.6485e4  # (coul)
        R = 8.3145 # (joule/degC)
        T = h.celsius + 273.19  # Kelvin
        E = (1e-3) * v  # convert mV to V
        Ci = ci + (self.soma().cap.monovalPerm) * (self.soma().cap.monovalConc)  #       : Monovalent permeability
        if (np.fabs(1-np.exp(-z*(F*E)/(R*T))) < 1e-6):  #denominator is small -> Taylor series
            ghk = (1e-6) * z * F * (Ci-co*np.exp(-z*(F*E)/(R*T)))*(1-(z*(F*E)/(R*T)))
        else:
            ghk = (1e-6) * z**2.*(E*F**2.)/(R*T)*(Ci-co*np.exp(-z*(F*E)/(R*T)))/(1-np.exp(-z*(F*E)/(R*T)))
        return ghk
