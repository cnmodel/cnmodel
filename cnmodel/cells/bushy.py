from __future__ import print_function
from neuron import h
from collections import OrderedDict
from .cell import Cell
from .. import synapses
from ..util import nstomho
from ..util import Params
import numpy as np
from .. import data
import pprint
pp = pprint.PrettyPrinter(indent=4, width=60)
    
__all__ = ['Bushy', 'BushyRothman']


class Bushy(Cell):
    
    celltype = 'bushy'
    spike_source = None
    scaled = False  # only allow scaling once
    
    @classmethod
    def create(cls, model='RM03', **kwds):
        if model == 'RM03':
            return BushyRothman(**kwds)
        else:
            raise ValueError ('Bushy model %s is unknown', model)

    def make_psd(self, terminal, psd_type, **kwds):
        """
        Connect a presynaptic terminal to one post section at the specified location, with the fraction
        of the "standard" conductance determined by gbar.
        The default condition is designed to pass the unit test (loc=0.5)
        
        Parameters
        ----------
        terminal : Presynaptic terminal (NEURON object)
        
        psd_type : either simple or multisite PSD for bushy cell
        
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
                self.NMDAR_vshift = data.get('sgc_synapse', species=self.species,
                        post_type=self.celltype, field='NMDAR_vshift')
                # adjust gmax to correct for initial Pr
                self.AMPAR_gmax = self.AMPAR_gmax/self.Pr
                self.NMDAR_gmax = self.NMDAR_gmax/self.Pr

#               original values (now in synapses.py):
#                self.AMPA_gmax = 3.314707700918133*1e3  # factor of 1e3 scales to pS (.mod mechanisms) from nS.
#                self.NMDA_gmax = 0.4531929783503451*1e3
                if 'AMPAScale' in kwds:  # normally, this should not be done!
                    self.AMPAR_gmax = self.AMPAR_gmax * kwds['AMPAScale']  # allow scaling of AMPA conductances
                if 'NMDAScale' in kwds:
                    self.NMDAR_gmax = self.NMDAR_gmax * kwds['NMDAScale']  # and NMDA... 
                return self.make_glu_psd(post_sec, terminal, self.AMPAR_gmax, 
                            self.NMDAR_gmax, loc=loc, nmda_vshift=self.NMDAR_vshift)
            elif terminal.cell.celltype == 'dstellate':
                return self.make_gly_psd(post_sec, terminal, psdtype='glyslow', loc=loc)
            elif terminal.cell.celltype == 'tuberculoventral':
                return self.make_gly_psd(post_sec, terminal, psdtype='glyslow', loc=loc)
            else:
                raise TypeError("Cannot make PSD for %s => %s" % 
                            (terminal.cell.celltype, self.celltype))
        else:
            raise ValueError("Unsupported psd type %s" % psd_type)

    def make_terminal(self, post_cell, term_type, **kwds):
        if term_type == 'simple':
            return synapses.SimpleTerminal(self.soma, post_cell, **kwds)

        elif term_type == 'multisite':
            if post_cell.celltype in ['mso']:
                nzones = data.get('bushy_synapse', species=self.species,
                        post_type=post_cell.celltype, field='n_rsites')
                delay = data.get('bushy_synapse', species=self.species,
                        post_type=post_cell.celltype, field='delay')
            else:
                raise NotImplementedError("No knowledge as to how to connect Bushy cell to cell type %s" %
                                        type(post_cell))
            pre_sec = self.soma
            return synapses.StochasticTerminal(pre_sec, post_cell, nzones=nzones,  spike_source=self.spike_source,
                                            delay=delay, **kwds)
        else:
            raise ValueError("Unsupported terminal type %s" % term_type)


class BushyRothman(Bushy):
    """
    VCN bushy cell models.
        Rothman and Manis, 2003abc (Type II, Type II-I)
        Xie and Manis, 2013
    """

    def __init__(self, morphology=None, decorator=None, nach=None,
                 ttx=False, species='guineapig', modelType=None, modelName=None,
                 debug=False, temperature=None):
        """
        Create a bushy cell, using the default parameters for guinea pig from
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
        
        species: string (default 'guineapig')
            species defines the pattern of ion channel densities that will be inserted, according to 
            prior measurements in various species. Note that
            if a decorator function is specified, this argument is ignored as the decorator will
            specify the channel density.

        modelName: string (default: None)
            modelName specifies the source conductance pattern (RM03, XM13, etc).
            modelName is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
                             
        modelType: string (default: None)
            modelType specifies the subtype of the cell model that will be used (e.g., "II", "II-I", etc).
            modelType is passed to the decorator, or to species_scaling to adjust point (single cylinder) models.
            
        debug: boolean (default: False)
            When True, there will be multiple printouts of progress and parameters.
            
        Returns
        -------
            Nothing
                 
        """
        super(BushyRothman, self).__init__()
        self.i_test_range={'pulse': (-1, 1, 0.05)}  # note that this might get reset with decorator according to channels
                                                    # The default values are set in the species_scaling routine
        if modelType == None:
            modelType = 'II'

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
            elif modelName  == 'XM13_nacncoop':
                dataset = 'XM13_nacncoop_channels'
            elif modelName  == 'XM13_nacn':
                dataset = 'XM13_nacn_channels'
            elif modelName  == 'XM13_nabu':
                dataset = 'XM13_nabu_channels'
            elif modelName.startswith('mGBC'):  # mouse Globular bushy - a bit different (and experimental)
                dataset = 'mGBC_channels'
            else:
                raise ValueError(f"ModelName {self.status['modelName']:s} not recognized for mouse {self.celltype:s} cells")
        else:
            raise ValueError(f"Species {species:s} not recognized for {self.celltype:s} cells")

        self.debug = debug
        self.status = {'species': species, 'cellClass': self.celltype, 'modelType': modelType, 'modelName': modelName,
                       'soma': True, 'axon': False, 'dendrites': False, 'pumps': False, 'hillock': False, 
                       'initialsegment': False, 'myelinatedaxon': False, 'unmyelinatedaxon': False,
                       'na': nach, 'ttx': ttx, 'name': self.celltype,
                       'morphology': morphology, 'decorator': decorator, 'temperature': temperature}
        self.vrange = [-70., -55.]  # set a default vrange for searching for rmp

        self.c_m = 0.9  # default in units of uF/cm^2
        self.spike_threshold = -40
        if self.debug:
            print( 'model type, model name, species: ', modelType, modelName, species, nach)


        self._valid_temperatures = (temp, )
        if self.status['temperature'] == None:
            self.status['temperature'] = temp

        soma = self.do_morphology(morphology)

        self.pars = self.get_cellpars(dataset, species=species, modelType=modelType)
        self.status['na'] = self.pars.natype
        # decorate the morphology with ion channels
        if decorator is None:   # basic "point" model, only on the soma, uses table data for soma.
            self.mechanisms = ['klt', 'kht', 'ihvcn', 'leak', self.pars.natype]
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ena = self.e_na
            self.soma.ek = self.e_k
            self.soma().ihvcn.eh = self.e_h
            self.soma().leak.erev = self.e_leak
            self.species_scaling(silent=True)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments, with tables.
            self.decorate()
        self.save_all_mechs()  # save all mechanisms inserted, location and gbar values...
        self.get_mechs(self.soma)
        if debug:
            print (f"   << Created {self.celltype.title():s} cell >>")

    def get_cellpars(self, dataset, species='guineapig', modelType='II'):
        """
        Retrieve parameters for the specifed model type and species from the data tables
        
        dataset : str (no default)
            name of the data table to use
        
        species : str (default: 'guineapig')
            Species table to use
        
        modelType : str (default: 'II')
            Model type to get parameters from the table.
        
        """
        cellcap = data.get(dataset, species=species, model_type=modelType,
            field='soma_Cap')  # total somatic capacitance (point cells)
        chtype = data.get(dataset, species=species, model_type=modelType,
            field='na_type')
        pars = Params(cap=cellcap, natype=chtype)
        if self.status['modelName'] == 'RM03':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'klt_gbar', 'ih_gbar', 'leak_gbar']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        elif self.status['modelName'] == 'XM13':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'klt_gbar', 'ihvcn_gbar', 'leak_gbar']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        elif self.status['modelName'] == 'XM13nacncoop':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'klt_gbar', 'ihvcn_gbar', 'leak_gbar']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        elif self.status['modelName'] == 'mGBC':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'klt_gbar', 'ihvcn_gbar', 'leak_gbar']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        elif self.status['modelName'] == 'XM13_nabu':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'klt_gbar', 'ihvcn_gbar', 'leak_gbar']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        elif self.status['modelName'] == 'XM13_nacncoop':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'klt_gbar', 'ihvcn_gbar', 'leak_gbar']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        elif self.status['modelName'] == 'XM13_nacn':
            for g in ['%s_gbar' % pars.natype, 'kht_gbar', 'klt_gbar', 'ihvcn_gbar', 'leak_gbar']:
                pars.additem(g,  data.get(dataset, species=species, model_type=modelType,
                    field=g))
        else:
            raise ValueError(f"get_cellpars: Model name {self.status['modelName']} is not yet implemented for cell type {self.celltype.title():s}")

        if self.debug:
            pars.show()
        return pars
        
    def species_scaling(self, silent=True):
        """
        This is called for POINT CELLS ONLY
        Adjust all of the conductances and the cell size according to the species requested.
        This scaling should be used ONLY for point models, as no other compartments
        are scaled.
        
        This scaling routine also sets the temperature for the model to a default value. Some models
        can be run at multiple temperatures, and so a default from one of the temperatures is used.
        The calling cell.set_temperature(newtemp) will change the conductances and reinitialize
        the cell to the new temperature settings.
        
        get_cellpars must be called before this is called.
        
        Parameters
        ----------
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        
        """
        assert self.scaled is False  # block double scaling!
        self.scaled = True
         
        soma = self.soma

        if self.status['species'] == 'mouse':
            # use conductance levels determined from Cao et al.,  J. Neurophys., 2007. as 
            # model description in Xie and Manis 2013. Note that
            # conductances were not scaled for temperature (rates were)
            # so here we reset the default Q10's for conductance (g) to 1.0
            if self.status['modelType'] not in ['II', 'II-I']:
                raise ValueError(f"\nModel type {self.status['modelType']:s} is not implemented for mouse {self.celltype.title():s} cells")
            if self.debug:
                print (f"  Setting conductances for mouse {self.celltype.title():s} cell ({self.status['modelType']})")

            self.vrange = [-68., -55.]  # set a default vrange for searching for rmp
            self.i_test_range = {'pulse': (-1., 1.0, 0.05)}
            self._valid_temperatures = (34., )
            if self.status['temperature'] is None:
                self.status['temperature'] = 34. 

            self.set_soma_size_from_Cm(self.pars.cap)
            self.adjust_na_chans(soma, sf=1.0)
            soma().kht.gbar = nstomho(self.pars.kht_gbar, self.somaarea)
            soma().klt.gbar = nstomho(self.pars.klt_gbar, self.somaarea)
            soma().ihvcn.gbar = nstomho(self.pars.ihvcn_gbar, self.somaarea)
            soma().leak.gbar = nstomho(self.pars.leak_gbar, self.somaarea)
            self.axonsf = 0.57
            
        elif self.status['species'] == 'guineapig':
            if self.debug:
                print ("  Setting conductances for guinea pig {self.celltype.title():s} {self.status['modelType']:s} cell, Rothman and Manis, 2003")

            self.i_test_range = {'pulse': (-0.4, 0.4, 0.02)}
            self.vrange = [-70., -55.]  # set a default vrange for searching for rmp
            self._valid_temperatures = (22., 38.)
            if self.status['temperature'] is None:
                self.status['temperature'] = 22. 

            sf = 1.0
            if self.status['temperature'] == 38.:  # adjust for 2003 model conductance levels at 38
                sf = 2  # Q10 of 2, 22->38C. (p3106, R&M2003c)
                # note that kinetics are scaled in the mod file.
            self.set_soma_size_from_Cm(self.pars.cap)
            self.adjust_na_chans(soma, sf=sf)
            soma().kht.gbar = nstomho(self.pars.kht_gbar, self.somaarea)
            soma().klt.gbar = nstomho(self.pars.klt_gbar, self.somaarea)
            soma().ihvcn.gbar = nstomho(self.pars.ih_gbar, self.somaarea)
            soma().leak.gbar = nstomho(self.pars.leak_gbar, self.somaarea)
            self.axonsf = 0.57
            
        else:
            errmsg = 'Species "%s" or model type "%s" is not recognized for Bushy cells.' %  (self.status['species'], self.status['modelType'])
            errmsg += '\n  Valid species are: \n'
            for s in knownspecies:
                errmsg += '    %s\n' % s
            errmsg += '-'*40
            raise ValueError(errmsg)

        self.check_temperature()
        if not silent:
           print (' set cell as: ', species)
           print (' with Vm rest = %6.3f' % self.vm0)

    def get_distancemap(self):
        return {'dend': {'klt': {'gradient': 'exp', 'gminf': 0., 'lambda': 50.},
                                 'kht': {'gradient': 'exp', 'gminf': 0., 'lambda': 50.},
                                 'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 50.}}, # linear with distance, gminf (factor) is multiplied by gbar
                        'dendrite': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                 'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                 'nav11': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                        'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                 'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                 'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 200.}}, # gradients are: flat, linear, exponential
                        }

    def add_axon(self):
        """
        Add a default axon from the generic cell class to the bushy cell (see cell class).
        """
        Cell.add_axon(self, self.c_m, self.R_a, self.axonsf)

    def add_pumps(self):
        """
        Insert mechanisms for potassium ion management, sodium ion management, and a
        sodium-potassium pump at the soma.
        """
        soma = self.soma
        soma.insert('k_conc')

        ki0_k_ion = 140
        soma().ki = ki0_k_ion
        soma().ki0_k_conc = ki0_k_ion
        soma().beta_k_conc = 0.075

        soma.insert('na_conc')
        nai0_na_ion = 5
        soma().nai = nai0_na_ion
        soma().nai0_na_conc = nai0_na_ion
        soma().beta_na_conc = 0.075

        soma.insert('nakpump')
        soma().nakpump.inakmax = 8
        soma().nao = 145
        soma().ko = 5
        soma().nakpump.Nai_inf = 5
        soma().nakpump.Ki_inf = 140
        soma().nakpump.ATPi = 5
        self.status['pumps'] = True

    def add_dendrites(self):
        """
        Add a simple dendrite to the bushy cell.
        """
        if self.debug:
            print ('Adding dendrite to Bushy model')
        section = h.Section
        primarydendrite = section(cell=self.soma)
        primarydendrite.connect(self.soma)
        primarydendrite.nseg = 10
        primarydendrite.L = 100.0
        primarydendrite.diam = 2.5
        primarydendrite.insert('klt')
        primarydendrite.insert('ihvcn')
        primarydendrite().klt.gbar = self.soma().klt.gbar / 2.0
        primarydendrite().ihvcn.gbar = self.soma().ihvcn.gbar / 2.0

        primarydendrite.cm = self.c_m
        primarydendrite.Ra = self.R_a
        nsecd = range(0, 5)
        secondarydendrite = []
        for ibd in nsecd:
            secondarydendrite.append(section(cell=self.soma))
        for ibd in nsecd:
            secondarydendrite[ibd].connect(primarydendrite)
            secondarydendrite[ibd].diam = 1.0
            secondarydendrite[ibd].L = 15.0
            secondarydendrite[ibd].cm = self.c_m
            secondarydendrite[ibd].Ra = self.R_a
        self.primarydendrite = primarydendrite
        self.secondarydendrite = secondarydendrite
        self.status['dendrite'] = True
        if self.debug:
            print ('Bushy: added dendrites')
            h.topology()
        self.add_section(maindend, 'primarydendrite')
        self.add_section(secdend, 'secondarydendrite')

