from neuron import h
import neuron as nrn
from ..util import nstomho
import numpy as np
from .cell import Cell
from .. import synapses
from .. import an_model
from .. import data

__all__ = ['SGC', 'SGC_TypeI', 'DummySGC']


class SGC(Cell):
    type = 'sgc'

    @classmethod
    def create(cls, model='I', species='mouse', **kwds):
        if model == 'dummy':
            return DummySGC(**kwds)
        elif model == 'I':
            return SGC_TypeI(species=species, **kwds)
        else:
            raise ValueError ('SGC model %s is unknown', model)
        
    def __init__(self, cf=None, sr=None):
        Cell.__init__(self)
        self._cf = cf
        self._sr = sr
        self.spike_source = None  # used by DummySGC to connect VecStim to terminal

    @property
    def cf(self):
        """ Center frequency
        """
        return self._cf
    
    @property
    def sr(self):
        """ Spontaneous rate group. 1=low, 2=mid, 3=high
        """
        return self._sr
        
    def make_terminal(self, post_cell, term_type, **kwds):
        """Create a StochasticTerminal and configure it according to the 
        postsynaptic cell type.
        """
        pre_sec = self.soma
        
        # Return a simple terminal unless a stochastic terminal was requested.
        if term_type == 'simple':
            return synapses.SimpleTerminal(pre_sec, post_cell, 
                                           spike_source=self.spike_source, **kwds)
        elif term_type == 'multisite':
            n_rsites = data.get('sgc_synapse', species='mouse', post_type=post_cell.type,
                            field='n_rsites')
            opts = {'nzones': n_rsites, 'delay': 0, 'dep_flag' : 1}
            opts.update(kwds)
            # when created, depflag is set True (1) so that we compute the DKR D*F to get release
            # this can be modified prior to the run by setting the terminal(s) so that dep_flag is 0
            # (no DKR: constant release probability)
            term = synapses.StochasticTerminal(pre_sec, post_cell,
                    spike_source=self.spike_source, **opts)
            
            kinetics = data.get('sgc_ampa_kinetics', species='mouse', post_type=post_cell.type,
                            field=['tau_g', 'amp_g'])
            term.set_params(**kinetics)
            dynamics = data.get('sgc_release_dynamics', species='mouse', post_type=post_cell.type,
                                field=['F', 'k0', 'kmax', 'kd', 'kf', 'taud', 'tauf', 'dD', 'dF'])
            term.set_params(**dynamics)
            return term
        else:
            raise ValueError("Unsupported terminal type %s" % term_type)
    

class DummySGC(SGC):
    """ SGC class with no cell body; this cell only replays a predetermined
    spike train.
    """
    def __init__(self, cf=None, sr=None, simulator='matlab'):
        """
        Parameters
        ----------
        cf : float (default: None)
            Required: the characteristic frequency for the SGC
        
        sr : int (default None)
            required : Selects the spontaneous rate group from the
            Zilany et al (2010) model. 1 = LSR, 2 = MSR, 3 = HSR
        
        simul(default: 'matlab')
            Sets the simulator interface that will be used. All models
            currently use the Zilany et al. model, but the simulator can
            be run though a Python-interface direcltly to the Matlab code
            as publicy available, (simulator='matlab'), or can be run through
            Rudieki's Python interface to the simulator's C code 
            (simulator='cochlea'). Requires installation of the modified
            versions of cochlea and thorns from github.com/pbmanis/cochlea and
            github.com/pbmanis/thorns.
        
        """
        SGC.__init__(self, cf, sr)
        self.vecstim = h.VecStim()
        
        # this causes the terminal to receive events from the VecStim:
        self.spike_source = self.vecstim
        
        # just an empty section for holding the terminal
        self.add_section(h.Section(), 'soma')
        
    def set_spiketrain(self, times):
        """ Set the times of spikes to be replayed by the cell.
        """
        self._spiketrain = times
        self._stvec = h.Vector(times)
        self.vecstim.play(self._stvec)

    def set_sound_stim(self, stim, seed, simulator='matlab'):
        """ Set the sound stimulus used to generate this cell's spike train.
        """
        self._sound_stim = stim
        spikes = an_model.get_spiketrain(cf=self.cf, sr=self.sr, seed=seed, 
            stim=stim, simulator=simulator)
        self.set_spiketrain(spikes * 1000)


class SGC_TypeI(SGC):
    """
    Spiral ganglion cell model
    
    """
    def __init__(self, morphology=None, decorator=None, nach=None, ttx=False,
                 species='guineapig', 
                 modelType='bm', cf=None, sr=None, debug=False):
        """
        Initialize a spiral ganglion Type I cell, based on the based on a bushy cell model.
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

        nach : string (default: 'na')
            nach selects the type of sodium channel that will be used in the model. A channel mechanim
            by that name must exist. The default is jsrna (Rothman et al., 1993)
        
        ttx : Boolean (default: False)
            If ttx is True, then the sodium channel conductance is set to 0 everywhere in the cell.
            Currently, this is not implemented.
        
        species: string (default 'guineapig')
            species defines the channel density that will be inserted for different models. Note that
            if a decorator function is specified, this argument is ignored.
            
        modelType: string (default: None)
            modelType specifies the type of the model that will be used. SGC model know about "a" (apical)
            and "bm" (basal-middle) models, based on Liu et al., JARO, 2014.
            modelType is passed to the decorator, or to species_scaling to adjust point models.

        cf : float (default: None)
            The CF for the auditory nerve fiber that this SGC represents.

        sr : string (default: None)
            The spontaneous rate group to which this fiber belongs. "LS", "MS", and "HS" are known values.

        debug: boolean (default: False)
            debug is a boolean flag. When set, there will be multiple printouts of progress and parameters.
            
        Returns
        -------
            Nothing
        
        """         
        
        super(SGC_TypeI, self).__init__(cf=cf, sr=sr)
        if modelType == None:
            modelType = 'bm'  # modelTypes are: a (apical), bm (basal middle)
        if nach == None:
            nach = 'jsrna'
        self.status = {'soma': True, 'axon': False, 'dendrites': False, 'pumps': False,
                       'na': nach, 'species': species, 'modelType': modelType, 'ttx': ttx, 'name': 'SGC',
                        'morphology': morphology, 'decorator': decorator, 'temperature': None}

        self.i_test_range={'pulse': [(-0.3, 0.3, 0.02), (-0.03, 0., 0.005)]}  # include finer range as well
        self.vrange = [-75., -55.]
        if morphology is None:
            """
            instantiate a basic soma-only ("point") model
            """
            soma = h.Section(name="SGC_Soma_%x" % id(self)) # one compartment of about 29000 um2
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
            self.mechanisms = [nach, 'klt', 'kht', 'leak']
            if modelType == 'a':
                self.mechanisms.append('ihsgcApical')
            elif modelType == 'bm':
                self.mechanisms.append('ihsgcBasalMiddle')
            else:
                raise ValueError ('Type %s not known for SGC model' % modelType)
            for mech in self.mechanisms:
                self.soma.insert(mech)
            self.soma.ek = self.e_k
            self.soma().leak.erev = self.e_leak
            self.species_scaling(silent=True, species=species, modelType=modelType)  # set the default type II cell parameters
        else:  # decorate according to a defined set of rules on all cell compartments
            self.decorate()
#        print 'Mechanisms inserted: ', self.mechanisms
        self.get_mechs(self.soma)
#        self.cell_initialize(vrange=self.vrange)
        if debug:
            print "<< SGC: Spiral Ganglion Cell created >>"

    def species_scaling(self, silent=True, species='guineapig', modelType='a'):
        """
        Adjust all of the conductances and the cell size according to the species requested.
        Used ONLY for point models.
        
        Parameters
        ----------
        species : string (default: 'guineapig')
            name of the species to use for scaling the conductances in the base point model
            Must be one of mouse or guineapig
        
        modelType: string (default: 'a')
            definition of HCN model type from Liu et al. JARO 2014:
            'a' for apical model
            'bm' for basal-middle model
        
        silent : boolean (default: True)
            run silently (True) or verbosely (False)
        
        Returns
        -------
        Nothing
        
        Notes
        -----
        The 'guineapig' model uses the mouse HCN channel model, verbatim. This may not
        be appropriate, given that the other conductances are scaled up.
        
        """
        
        soma = self.soma
        if species == 'mouse':
            self.set_soma_size_from_Cm(12.0)
            self._valid_temperatures = (34.,)
            if self.status['temperature'] is None:
                self.set_temperature(34.)
            self.adjust_na_chans(soma, gbar=350.)
            soma().kht.gbar = nstomho(58.0, self.somaarea)
            soma().klt.gbar = nstomho(80.0, self.somaarea)
                # nstomho(200.0, somaarea) * scalefactor
            if modelType == 'a':
                soma().ihsgcApical.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcApical.eh = -41
            elif modelType == 'bm':
                soma().ihsgcBasalMiddle.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcBasalMiddle.eh = -41
            else:
                raise ValueError('Ihsgc modelType %s not recognized for species %s' % (species, modelType))
            soma().leak.gbar = nstomho(2.0, self.somaarea)

        elif species == 'guineapig':
            # guinea pig data from Rothman and Manis, 2003, modelType II
            self.set_soma_size_from_Cm(12.0)
            self._valid_temperatures = (22.,)
            if self.status['temperature'] is None:
                self.set_temperature(22.)
            self.adjust_na_chans(soma, gbar=1000.)
            soma().kht.gbar = nstomho(150.0, self.somaarea)
            soma().klt.gbar = nstomho(200.0, self.somaarea)
                # nstomho(200.0, somaarea) * scalefactor
            if modelType == 'a':
                soma().ihsgcApical.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcApical.eh = -41
            elif modelType == 'bm':
                soma().ihsgcBasalMiddle.gbar = nstomho(3.0, self.somaarea)
                soma().ihsgcBasalMiddle.eh = -41
            else:
                raise ValueError('Ihsgc modelType %s not recognized for species %s' % (species, modelType))
            soma().leak.gbar = nstomho(2.0, self.somaarea)

        else:
            raise ValueError('Species %s or species-modelType %s is not recognized for SGC cells' % (species, modelType))

        self.status['species'] = species
        self.status['modelType'] = modelType
#        self.cell_initialize(showinfo=False)
        self.check_temperature()
        if not silent:
            print 'set cell as: ', species
            print ' with Vm rest = %f' % self.vm0

    def adjust_na_chans(self, soma, gbar = 1000., debug=False):
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
                print "sgc using inva11"
            print 'nav11 gbar: ', soma().nav11.gbar
        elif nach in ['na', 'nacn']:
            soma().na.gbar = gnabar
            soma.ena = self.e_na
            if debug:
                print 'na gbar: ', soma().na.gbar
        else:
            raise ValueError('Sodium channel %s is not recognized for SGC cells', nach)

    def i_currents(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point
        vrange brackets the interval
        Implemented here are the basic RM03 mechanisms
        This function should be replaced for specific cell types.
        """
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = V

        h.t = 0.
        h.celsius = self.status['temperature']
        h.finitialize()
        self.ix = {}
        if 'na' in self.mechanisms:
            #print dir(self.soma().na)
            self.ix['na'] = self.soma().na.gna*(V - self.soma().ena)
        if 'jsrna' in self.mechanisms:
            #print dir(self.soma().na)
            self.ix['jsrna'] = self.soma().jsrna.gna*(V - self.soma().ena)
        if 'klt' in self.mechanisms:
            self.ix['klt'] = self.soma().klt.gklt*(V - self.soma().ek)
        if 'kht' in self.mechanisms:
            self.ix['kht'] = self.soma().kht.gkht*(V - self.soma().ek)
        if 'ihsgcApical' in self.mechanisms:
            self.ix['ihsgcApical'] = self.soma().ihsgcApical.gh*(V - self.soma().ihsgcApical.eh)
        if 'ihsgcBasalMiddle' in self.mechanisms:
            self.ix['ihsgcBasalMiddle'] = self.soma().ihsgcBasalMiddle.gh*(V - self.soma().ihsgcBasalMiddle.eh)
        if 'leak' in self.mechanisms:
            self.ix['leak'] = self.soma().leak.gbar*(V - self.soma().leak.erev)
#        print self.status['name'], self.status['type'], V, self.ix
        return np.sum([self.ix[i] for i in self.ix])

