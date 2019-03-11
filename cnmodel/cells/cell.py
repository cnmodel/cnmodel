from __future__ import print_function
import weakref
import numpy as np
import scipy.optimize
from collections import OrderedDict
import neuron
from neuron import h
from ..util import nstomho, mho2ns
from ..util import custom_init
from .. import synapses
from .. import data
from .. import morphology
from .. import decorator

"""
Term definitions:
cell class is the class of morphological cell: bushy, tstellate, etc. 
Each cell class is implmeneted as a separate python class (no pun)
modelName is name of the source model used so it is like the type, but one level up). 
ModelNames are RM03, XM13, and for other cell types may refer to the original model, 
such as POK (Kanold pyramidal cell), MCG (McGinley octopus), Eager, etc.
These model designations may have only one model type (POK), or may have multiple types (RM03, XM13)
modelType refers to the Rothman and Manis 2003 model classes (I, II, I-c, I-t, II-1, I-II, etc)
These are physiologically based, but in the ion channel tables are mapped to morphological classes sort of,


"""

class Cell(object):
    """
    Base class for all cell types.
    
    """
    type = None
    
    # create a lookup table to map sections to their parent cell
    sec_lookup = weakref.WeakValueDictionary()
    
    @classmethod
    def from_section(cls, sec):
        return cls.sec_lookup[sec.name()]
    
    def __init__(self):
        # dictionary of all sections associated with this cell
        self.hr = None # hoc reader - e.g., we have read a morphology file.
        self.all_sections = {}
        # the following section types (parts) are known to us:
        for k in ['soma', 'maindend', 'secdend', 'dend', 'dendrite', 'primarydendrite', 'secondarydendrite',
            'internode',
            'initialsegment', 'axonnode', 'axon', 'unmyelinatedaxon', 'myelinatedaxon', 'hillock']:
            self.all_sections[k] = []  # initialize to an empty list
        self.species = 'mouse'
        self.status = {}  # dictionary of parameters used to instantiate the cell.
        # Record synaptic inputs and projections
        self.inputs = []  # inputs are recorded - synapse object, post_opts and kwds
        self.outputs = []
        self.initial_mechanisms = None
        # each cell has the following parameters:
        self.totcap = None  # total membrane capacitance (somatic)
        self.somaarea = None  # total soma area
        self.initsegment = None  # hold initial segment sections
        self.axnode = None  # hold nodes of ranvier sections
        self.internode = None  # hold internode sections
        self.maindend = None  # hold main dendrite sections
        self.secdend = None  # hold secondary dendrite sections
        self.dendrite = None
        self.axon = None
        self.axonsf = None  # axon diameter scale factor
        # define defaults for these parameters (RM03 model defaults)
        self.e_k = -70  # potassium reversal potential, mV
        self.e_na = 55
        self.e_h = -43
        self.c_m = 0.9  # specific membrane capacitance,  uf/cm^2
        self.R_a = 150  # axial resistivity of cytoplasm/axoplasm, ohm.cm
        self.e_leak = -65
        # Recommended current (min, max, step) for testing this cell
        self.i_test_range=(-0.5, 0.5, 0.05)  # defines default current steps for IC curve
        
        # Recommended threshold for detecting spikes from this cell
        self.spike_threshold = -40
        
        # Resting potential for this cell, determined by calling
        # self.find_i0()
        self.vm0 = None

    def check_temperature(self):
        if self.status['temperature'] not in self._valid_temperatures:
            tstring = ', '.join('%3.1f ' % t for t in self._valid_temperatures)
            raise ValueError('Cell %s %s %s temperature %3.1f is invalid; must be in: [%s]' %
                 (self.type, self.status['species'], self.status['modelType'], self.status['temperature'], tstring))

    def set_temperature(self, temperature):
        """
        Set the temperature setting for this cell. 
        """
        if self.status['decorator'] is None:
            if self.status['temperature'] is None:  # only if not already set
                self.status['temperature'] = temperature
                self.species_scaling(species=self.status['species'], modelType=self.status['modelType'])
        else:
            self.status['temperature'] = temperature
#            self.decorate()  # call the decorator
            
    def set_morphology(self, morphology_file=None):
        """
        Set the cell's morphological structure from a file that defines sections
        (for example, a morphology file read by neuronvis), or from a morphology
        object that has already been retrieved/created.
        
        Parameters
        ----------
        morphology_file : string or morphology object (default: None)
            File name/path for the morphology file (for example, .hoc or .swc file)
            Alternatively, this can be a morphology object returned by the morphology class.

        Returns
        -------
        nothing
            
        """
        self.morphology_file = morphology_file  # save the source file name
        if isinstance(morphology_file, str):
            if morphology_file.endswith('.hoc'):
                self.morphology = morphology.HocReader(morphology_file)
            elif morphology_file.endswith('.swc'):
                self.morphology = morphology.SwcReader(morphology_file)
            else:
                raise ValueError('Unknown morphology file type [must be .hoc or .swc]')
        elif isinstance(morphology_file, morphology.Morphology):
            self.morphology = morphology_file
        else:
            print(morphology_file)
            raise TypeError('Invalid morphology type: must be filename(str) or ')
        self.hr = self.morphology # extensive renaming required in calling classes, temporary fix.
        self.morphology.read_section_info()  # not sure this is necessary... 
        # these were not instantiated when the file was read, but when the decorator was run.
        for s in self.hr.sec_groups.keys():
            for sec in self.hr.sec_groups[s]:
                section =self.hr.get_section(sec)
                mechs = self.hr.get_mechanisms(sec)
                if s == 'myelinatedaxon':
                    section.cm = 0.002
                self.add_section(section, s) # add the section to the cell.
               # print '\nmechanisms for section: %s', section
               # self.print_mechs(section)
        self.set_soma_size_from_Section(self.soma)  # this is used for reporting and setting g values...
        if isinstance(self.soma, list):
            self.distances(self.soma[1])
        else:
            self.distances(self.soma)
        self.hr.distanceMap = self.distanceMap

    def add_section(self, sec, sec_type):
        """
        Add a section (or list of sections) to this cell. 
        This adds the section to self.all_sections[sec_type] and also allows 
        the cell to be accessed from the section using 
        cells.cell_from_section().
        
        Notes:
        
        *sec_type* must be one of the keys already in self.all_sections.
        
        This method does not connect sections together; that must be 
        done manually.
        
        """
        if not isinstance(sec, list):
            sec = [sec]
        self.all_sections[sec_type].extend(sec)
        for s in sec:
            Cell.sec_lookup[s.name()] = self
    
    def list_sections(self):
        # print self.all_sections
        print('Known Section names:')
        for sec in self.all_sections:
            print('  %s' % sec)
            s = self.all_sections[sec]
            # print 's: ', s
            if len(s) > 0:
                print('    ------------------------------------------')
                print('    Sections present:')
                for u in s:
                    print('    Type: %s (%s, %s): %s' % (sec,  u.name(), str(self.hr.get_section(u.name())), Cell.sec_lookup[u.name()]))
                print('    ------------------------------------------')
            else:
                print('    No section of this type in cell')
    
    def get_section_type(self, sec):
        for s in self.all_sections:
            if sec in self.all_sections[s]:
                return s
        return None

    def get_post_sec(self, kwds):
        """
        Get the postsynaptic section from the value of postsite 
        in kwds. This is typically called from the cell-specific make_psd method.
        If the key 'postsite' is in the kwds dict, we look it up. 
        If not, then we use the soma section as a default instead.
        
        Parameters
        ----------
        kwds : dict
             dictionary of keywords, may have a key 'postsite'
        
        Returns:
            loc, post_sec
                the location (0-1) of the desired point process insertion, and
                post_sec, the neuron section where that insertion will take place
        """
        if 'postsite' in kwds:  # use a defined location instead of the default (soma(0.5)
            postsite = kwds['postsite']
            loc = postsite[1]  # where on the section?
            uname = 'sections[%d]' % postsite[0]  # make a name to look up the neuron section object
            post_sec = self.hr.get_section(uname)  # Tell us where to put the synapse.
        else:
            loc = 0.5
            post_sec = self.soma
        return loc, post_sec

    def set_d_lambda(self, freq=100, d_lambda=0.1):
        """
        Sets nseg in each section to an odd value so that its segments are no longer than
        d_lambda x the AC length constant at frequency freq in that section.
        The defaults are reasonable values for most models
        Be sure to specify your own Ra and cm before calling geom_nseg()

        To understand why this works,
        and the advantages of using an odd value for nseg,
        see  Hines, M.L. and Carnevale, N.T. NEURON: a tool for neuroscientists. The Neuroscientist 7:123-135, 2001.
        This is a python version of the hoc code.
        
        Parameters
        ----------
        freq : float, default=100. (Hz)
            Frequency in Hz to use in computing nseg.
        d_lambda : float, default=0.1
            fraction of AC length constant for minimum segment length
        
        """
        if self.hr is None:  # no hoc reader file, so no adjustments
            return
        for st in self.all_sections.keys():
            for i, section in enumerate(self.all_sections[st]):
                nseg  = int((section.L/(d_lambda*self._lambda_f(freq, section))+0.9)/2)*2 + 1
                if nseg < 3:
                    nseg = 3 # ensure at least 3 segments per section...
                section.nseg = nseg


    def _lambda_f(self, freq, section):
        """
        get lambda_f for the section (internal)
        
        Parameters
        ----------
        freq : float, default=100. (Hz)
            Frequency in Hz to use in computing nseg.
        section : Neuron section object
        
        Returns
        -------
        section length normalized by the length constant at freq.
        """
        self.hr.h('access %s' % section.name())
        if self.hr.h.n3d() < 2:
            return 1e-5*np.sqrt(section.diam/(4.0*np.pi*freq*section.Ra*section.cm))
        # above was too inaccurate with large variation in 3d diameter
        # so now we use all 3-d points to get a better approximate lambda
        x1 = self.hr.h.arc3d(0)
        d1 = self.hr.h.diam3d(0)
        lam = 0.001
        for i in range(int(self.hr.h.n3d())-1):
                x2 = self.hr.h.arc3d(i)
                d2 = self.hr.h.diam3d(i)
                lam = lam + ((x2 - x1)/np.sqrt(d1 + d2))
                x1 = x2
                d1 = d2
        
        #  length of the section in units of lambda
        lam = lam * np.sqrt(2.0) * 1e-5*np.sqrt(4.0*np.pi*freq*section.Ra*section.cm)
        return section.L/lam


    @property
    def soma(self):
        """
        First (or only) section in the "soma" section group.
        """
        if isinstance(self.all_sections['soma'], list):
            return self.all_sections['soma'][0]
        else:
            return self.all_sections['soma']

    def decorate(self):
        """
        decorate the cell with it's own class channel decorator
        """
        self.decorated = decorator.Decorator(cell=self)
        self.decorated.channelValidate(self, verify=False)
        self.mechanisms = self.hr.mechanisms  # copy out all of the mechanisms that were inserted

    # def channel_manager(self, modelType='RM03'):
    #     """
    #     Every cell class should have a channel manager if it is set up to handle morphology.
    #     This function should be overridden in the class with an appropriate routine that
    #     builds the dictionary needed to decorate the cell. See the bushy cell class for
    #     an example.
    #
    #     Parameters
    #     ----------
    #     modelType : string (default: 'RM03')
    #          A string that identifies what type of model the channel manager will implement.
    #          This may be used to define different kinds of channels, or channel densities
    #          and compartmental placement for different cells.
    #     """
    #     raise NotImplementedError("No channel manager exists for cells of the class: %s" %
    #                               (self.__class__.__name__))

    def connect(self, post_cell, pre_opts=None, post_opts=None, **kwds):
        """
        Create a new synapse connecting this cell to a postsynaptic cell. 
        The synapse is automatically created using 
        pre_cell.make_terminal(post_cell, \**pre_opts) and  
        post_cell.make_psd(terminal, \**post_opts).
        
        By default, the cells decide which sections to connect. This can be 
        overridden by specifying 'section' in pre_opts and/or post_opts.
       
        Parameters
        ----------
        post_cell : NEURON section (required)
            The postsynaptic cell that will receive the connection.
        pre_opts : dictionary of options for the presynaptic cell (default: None)
            see the synapses class for valid options and format.
        post_opts : diction of options for the postsynaptic cell (default: None)
            see synapses class for valid options and format.
        \**kwds : (optional)
            argmuments that are passed to the synapses class.
        
        Returns
        -------
        the synapse object
        
        """
        if pre_opts is None:
            pre_opts = {}
        if post_opts is None:
            post_opts = {}
        
        synapse = synapses.Synapse(self, pre_opts, post_cell, post_opts, **kwds)
        self.outputs.append(synapse)
        post_cell.inputs.append([synapse, post_opts, kwds])
        
        return synapse

    def print_connections(self):
        """
        This is mostly for debugging ...
        """
        print('outputs: ', self.outputs)
        print('inputs: ', self.inputs)
    
    def make_terminal(self, post_cell, **kwds):
        """
        Create a synaptic terminal release mechanism suitable for output
        from this cell to post_sec
        This routine is a placeholder and should be replaced in the specific
        cell class with code that performs the required actions for that class.
        
        Parameters
        ----------
        post_cell : the target terminal cell (required)
        
        \**kwds : parameters passed to the terminal
        
        """
        raise NotImplementedError("Cannot make Terminal connecting %s => %s" % 
                                  (self.__class__.__name__, 
                                   post_cell.__class__.__name__))

    def make_psd(self, terminal, **kwds):
        """
        Create a PSD suitable for synaptic input from pre_sec.
        This routine is a placeholder and should be overridden in the specific
        cell class with code that performs the required actions for that class.
        
        Parameters
        ----------
        terminal : the terminal that connects to the PSD (required)
        
        \**kwds : parameters passed to the terminal
        
        """
        pre_cell = terminal.cell
        raise NotImplementedError("Cannot make PSD connecting %s => %s" %
                                  (pre_cell.__class__.__name__, 
                                   self.__class__.__name__))

    def make_glu_psd(self, post_sec, terminal, AMPA_gmax, NMDA_gmax, **kwds):
        # Get AMPAR kinetic constants from database 
        params = data.get('sgc_ampa_kinetics', species=self.species, post_type=self.type,
                            field=['Ro1', 'Ro2', 'Rc1', 'Rc2', 'PA'])
        
        return synapses.GluPSD(post_sec, terminal,
                                ampa_gmax=AMPA_gmax,
                                nmda_gmax=NMDA_gmax,
                                ampa_params=dict(
                                    Ro1=params['Ro1'],
                                    Ro2=params['Ro2'],
                                    Rc1=params['Rc1'],
                                    Rc2=params['Rc2'],
                                    PA=params['PA']),
                                **kwds)

    def make_gly_psd(self, post_sec, terminal, type, **kwds):
        # Get GLY kinetic constants from database 
        params = data.get('gly_kinetics', species=self.species, post_type=self.type,
                            field=['KU', 'KV', 'XMax'])
        psd = synapses.GlyPSD(post_sec, terminal,
                                psdType=type,
                                **kwds)
        return psd

    def make_exp2_psd(self, post_sec, terminal, weight=0.01, loc=0.5):
        return synapses.Exp2PSD(post_sec, terminal, weight=weight, loc=loc)

    def print_status(self):
        print("\nCell model: %s" % self.__class__.__name__)
        print(self.__doc__)
        print('    Model Status:')
        print('-'*24)
        for s in self.status.keys():
            print('{0:>12s} : {1:<12s}'.format(s, repr(self.status[s])))
        print('-'*32)

    def cell_initialize(self, showinfo=False, vrange=None, **kwargs):
        """
        Initialize this cell to it's "rmp" under current conditions
        All sections in the cell are set to the same value
        """
        if self.vm0 is None:
            self.vm0 = self.find_i0(showinfo=showinfo, vrange=vrange, **kwargs)
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = self.vm0

    def get_mechs(self, section):
        """
        return a list of the mechanisms that are present in a section
        a mechanism is required to have a gbar variable.
        This routine should be called at the end of every cell creation routine.
        """
        u = dir(section())
        mechs = []
        for m in u:
            if m[0:2] == '__':
                continue
            if m in ['cm', 'diam', 'k_ion', 'na_ion', 'next', 'point_processes', 'sec', 'v', 'x']:
                continue  # skip non-mechanisms known to us
            try:
                gx=eval('section().'+m+'.gbar')
                mechs.append(m)
            except:
                pass
        self.mechs = mechs
        return mechs

    def print_mechs(self, section):
        """
        print the mechanisms that are inserted into the specified section,
        and their densities (in uS/cm^2)
        
        """
        print('\n    Installed mechanisms:')
        self.get_mechs(section)
        #print eval('section().nav11.gbar')
 
        print('somaarea: {:.3e}'.format(self.somaarea))
        print('Mechanisms:', end='')
        for s in self.mechs:
            print(' {:>8s} '.format(s), end='')
        print('')
        for m in self.mechs:
            try:
                gx=eval('section().'+m+'.gbar')
                erev = 0.
                if m == 'leak':
                    erev = eval('section().'+m+'.erev')
                if m in ['jsrna', 'na', 'nacn', 'nav11', 'nacncoop', 'napyr', 'nap']:
                    erev = eval('section().ena')
                if m in ['klt', 'kht', 'ka']:
                    erev = eval('section().e%s'%m)
                if m in ['kis', 'kif', 'kdpyr', 'kcnq']:
                    erev = eval('section().ek')
                if m in ['hcno', 'ihvcn', 'hcnobo', 'ihpyr', 'ihpyr_adj']:
                    erev = eval('section().'+m+'.eh')
                print('{0:>12s} : {2:8.1f} nS  {1:7.3e} mho/cm2  {3:>5.1f} mV'.
                        format(m, gx, mho2ns(gx, self.somaarea), erev))
            except:
                print('{0:>12s} : <no gbar> '.format(m))
        print('-'*32)
        
    def print_all_mechs(self):
        print(self.get_all_mechs())

    def get_all_mechs(self):
        """
        return a string with all the mechanisms
        """
        res = '\nAll mechanisms in all sections: \n'
        for part in self.all_sections.keys():
            if len(self.all_sections[part]) == 0:
                #res += 'Cell part: %s hs not sections' % part
                continue
            res += 'Cell part: %s\n' % part
            for sec in self.all_sections[part]:
                res += '   Section: %s\n' % sec.name()
                res += '        %s' % self.get_mechs(sec) + '\n'
                for m in self.get_mechs(sec):
                    gx = eval('sec().'+m+'.gbar')
                    res += '            %s: %f\n' % (m, gx)
        return res

    def save_all_mechs(self):
        """
        get and save all of the initial mechanisms and their
        maximal conductances when the cell is created.
        We use this to get and check values later when the run
        is actually done.
        Note: some cell constructions may require that save_all_mechs
        be done again after the initial "build". In this case,
        setting the cell's initial_mechanisms property to None must
        be done to allow a new configuration of mechanisms to be saved.
        
        """
        if self.initial_mechanisms is not None:
            raise ValueError('Cells: Attempting to save initial mechanisms more than once')
        self.initial_mechanisms = {}
        for part in self.all_sections.keys():
            self.initial_mechanisms[part] = {}
#            print('Cell part: %s' % part )
            for sec in self.all_sections[part]:
#                print('   Section: ', sec)
#                print('        ', self.get_mechs(sec))
                self.initial_mechanisms[part][sec] = {}
                for m in self.get_mechs(sec):
                    gx = eval('sec().'+m+'.gbar')
#                    print('            %s: %f' % (m, gx))
                    self.initial_mechanisms[part][sec][m] = gx

    def check_all_mechs(self):
        """
        Check that all mechanisms are the same as when we initially created the cell
        """
        check = {}
        for part in self.all_sections.keys():
            if part not in self.initial_mechanisms.keys():
                raise ValueError('Cell part %s was not in the original cell')
            check[part] = {}
            for sec in self.all_sections[part]:
                #print('   Section: ', sec)
                #print('        ', self.get_mechs(sec))
                if sec not in self.initial_mechanisms[part].keys():
                    raise ValueError('Cell section was not in the original cell: ', sec)
                check[part][sec] = sec
                for m in self.get_mechs(sec):
                    gx = eval('sec().'+m+'.gbar')
                    #print('            %s: %f' % (m, gx))
                    if m not in self.initial_mechanisms[part][sec].keys():
                        raise ValueError('Mechanism %s was not in cell part %s, section = ' % (m, part), sec)
                    if self.initial_mechanisms[part][sec][m] != gx:
                        raise ValueError('Conductance for mechanism %s in cell part %s has changed (%f, %f), section = ' %
                            (m, part, self.initial_mechanisms[part][sec][m], gx), sec)
        return True
        
    def get_cellpars(self, dataset, species='guineapig', cell_type='II'):
        raise NotImplementedError('get_cellpars should be reimplemented in the individual cell modules')
    
    def channel_manager(self, modelName=None, modelType=None):
        """
        This routine defines channel density maps and distance map patterns
        for each type of compartment in the cell. The maps
        are used by the ChannelDecorator class (specifically, its private
        \_biophys function) to decorate the cell membrane.
        These settings are only used if the decorator is called; otherwise
        for point cells, the species_scaling routine defines the channel
        densities.
        
        Parameters
        ----------
        modelType : string (default: 'None'
            A string that defines the type of the model. 
            These are determined in the tables in the data directory, for ionchannels.py

        Returns
        -------
        Nothing
        
        Notes
        -----
        This routine defines the following variables for the class:
        
            * conductances (gBar)
            * a channelMap (dictonary of channel densities in defined anatomical compartments)
            * a current injection range for IV's (used for testing)
            * a distance map, which defines how each conductance in a selected compartment
              changes with distance from the soma. The current implementation includes both
              linear and exponential gradients,
              the minimum conductance at the end of the gradient, and the space constant or
              slope for the gradient.
        
        """

        dataset = '%s_channels' % modelName
        decorationmap = dataset + '_compartments'
        # print('dataset: {0:s}   decorationmap: {1:s}'.format(dataset, decorationmap))
        cellpars = self.get_cellpars(dataset, species=self.status['species'], modelType=modelType)
        refarea = 1e-3*cellpars.cap / self.c_m
# ?       print ('cellpars: ' )
        cellpars.show()
        # print(' species: ', self.status['species'])
        # print('m# odelType: ', modelType)
#         print('dataset: ', dataset)
        table = data._db.get_table_info(dataset)
#         table = data.get_table_info('mGBC_channels')
#         print(dir(data.ionchannels))
#         print( data.print_table('mGBC_channels'))
        
        if len(table.keys()) == 0:
            raise ValueError('data table %s lacks keys - does it exist?' % dataset)
        chscale = data._db.get_table_info(decorationmap)
        pars = {}
        # retrive the conductances from the data set
        # print ('table keys: ', table.keys())
        # print('table: ', table)
        # print('chscale: ', chscale)
        for g in table['field']:
            x = data._db.get(dataset, species=self.status['species'], model_type=modelType,
                                field=g)
            if not isinstance(x, float):
                continue
            if '_gbar' in g:
                pars[g] = x/refarea
            else:
                pars[g] = x
        
        self.channelMap = OrderedDict()
        for c in chscale['compartment']:
            self.channelMap[c] = {}
            for g in pars.keys():
                if g not in chscale['parameter']:
#                    print ('Parameter %s not found in chscale parameters!' % g)
                    continue
                scale = data._db.get(decorationmap, species=self.status['species'], model_type=modelType,
                        compartment=c, parameter=g)
                if '_gbar' in g:
                    self.channelMap[c][g] = pars[g]*scale
                else:
                    self.channelMap[c][g] = pars[g]

        self.irange = np.linspace(-0.6, 1, 9)
        self.distMap =         {'dend': {'klt': {'gradient': 'exp', 'gminf': 0., 'lambda': 50.},
                                 'kht': {'gradient': 'exp', 'gminf': 0., 'lambda': 50.},
                                 'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 50.}}, # linear with distance, gminf (factor) is multiplied by gbar
                        'dendrite': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                 'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                 'nav11': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.}}, # linear with distance, gminf (factor) is multiplied by gbar
                        'apic': {'klt': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                 'kht': {'gradient': 'linear', 'gminf': 0., 'lambda': 100.},
                                 'nav11': {'gradient': 'exp', 'gminf': 0., 'lambda': 200.}}, # gradients are: flat, linear, exponential
                        }
        self.check_temperature()
        return        




    def i_currents(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point
        vrange brackets the interval
        Implemented here are the basic known mechanisms. If you add or need
        more mechanisms, they either need to be accomadated in this routine,
        or this routine needs to be implemented (overridden) in the
        specific cell class. 
        
        """
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = V
        h.celsius = self.status['temperature']
        h.t = 0.
        h.finitialize(V)
        h.fcurrent()
        self.ix = {}

        if 'na' in self.mechanisms:
            #print dir(self.soma().na)
            try:
                self.ix['na'] = self.soma().na.gna*(V - self.soma().ena)
            except:
                self.ix['na'] = self.soma().nav11.gna*(V - self.soma().ena)
        if 'jsrna' in self.mechanisms:
            self.ix['jsrna'] = self.soma().jsrna.gna*(V - self.soma().ena)
        if 'nav11' in self.mechanisms:
            self.ix['nav11'] = self.soma().nav11.gna*(V - self.soma().ena)
        if 'nacn' in self.mechanisms:
            self.ix['nacn'] = self.soma().nacn.gna*(V - self.soma().ena)
        if 'napyr' in self.mechanisms:
            self.ix['napyr'] = self.soma().napyr.gna*(V - self.soma().ena)
        if 'nap' in self.mechanisms:
            self.ix['nap'] = self.soma().nap.gna*(V - self.soma().ena)
        if 'nacncoop' in self.mechanisms:
            self.ix['nacncoop'] = self.soma().nacncoop.gna*(V - self.soma().ena)
        
        if 'klt' in self.mechanisms:
            self.ix['klt'] = self.soma().klt.gklt*(V - self.soma().ek)
        if 'kht' in self.mechanisms:
            self.ix['kht'] = self.soma().kht.gkht*(V - self.soma().ek)
        if 'ka' in self.mechanisms:
            self.ix['ka'] = self.soma().ka.gka*(V - self.soma().ek)
        if 'kdpyr' in self.mechanisms:
            self.ix['kdpyr'] = self.soma().kdpyr.gk*(V - self.soma().ek)
        if 'kcnq' in self.mechanisms:
            self.ix['kcnq'] = self.soma().kdcnq.gk*(V - self.soma().ek)
        if 'kis' in self.mechanisms:
            self.ix['kis'] = self.soma().kis.gk*(V - self.soma().ek)
        if 'kif' in self.mechanisms:
            self.ix['kif'] = self.soma().kif.gk*(V - self.soma().ek)
        
        if 'ihvcn' in self.mechanisms:
            self.ix['ihvcn'] = self.soma().ihvcn.gh*(V - self.soma().ihvcn.eh)
        if 'ihpyr' in self.mechanisms:
            self.ix['ihpyr'] = self.soma().ihpyr.gh*(V - self.soma().ihpyr.eh)
        if 'ihpyr_adj' in self.mechanisms:
            self.ix['ihpyr_adj'] = self.soma().ihpyr_adj.gh*(V - self.soma().ihpyr_adj.eh)
        if 'hcno' in self.mechanisms:
            raise ValueError('HCNO is not supported - use hcnobo instead')
            #self.ix['hcno'] = self.soma().hcno.gh*(V - self.soma().hcno.eh)
        if 'hcnobo' in self.mechanisms:
            self.ix['hcnobo'] = self.soma().hcnobo.gh*(V - self.soma().hcnobo.eh)
        
        if 'leak' in self.mechanisms:
            self.ix['leak'] = self.soma().leak.gbar*(V - self.soma().leak.erev)
#        print self.status['name'], self.status['type'], V, self.ix
        isum = np.sum([self.ix[i] for i in self.ix])
#        print 'conductances: ', self.ix.keys()
#        print 'V, isum, values: ', V, isum, [self.ix[i] for i in self.ix]
        return isum

    def find_i0(self, vrange=None, showinfo=False):
        """
        find the root of the system of equations in vrange.
        Finds RMP fairly accurately as zero current level for current conductances.
        
        Parameters
        ----------
        vrange : list of 2 floats (default: [-70, -55])
            The voltage range over which the root search will be performed.
            
        showinfo : boolean (default: False)
            a flag to print out which roots were found and which mechanisms were in the cell
            
        Returns
        -------
        The voltage at which I = 0 in the vrange specified
        """
        if vrange is None:
             vrange = self.vrange
        # print( vrange)
        # print (self.i_currents(V=vrange[0]), self.i_currents(V=vrange[1]))
        # v0 = scipy.optimize.brentq(self.i_currents, vrange[0], vrange[1], maxiter=10000)
        # print( 'v0: ', v0)
        try:
            v0 = scipy.optimize.brentq(self.i_currents, vrange[0], vrange[1], maxiter=10000)
        except:
            print('find i0 failed:')
            print(self.ix)
            i0 = self.i_currents(V=vrange[0])
            i1 = self.i_currents(V=vrange[1])
            ivi = []
            ivv = []
            for v in np.arange(vrange[0], vrange[1], 0.5):
                ivi.append(self.i_currents(V=v))
                ivv.append(v)
            print ('iv: ')
            for i in range(len(ivi)):
                print('%6.1f  %9.4f' % (ivv[i], ivi[i]))
            print('This means the voltage range for the search might be too large\nor too far away from the target')
            raise ValueError('vrange not good for %s : %f at %6.1f, %f at %6.1f, temp=%6.1f' %
                             (self.status['name'], i0, vrange[0], i1, vrange[1], h.celsius))
        # check to be sure all the currents that are needed are calculated
        # can't do this until i_currents has populated self.ix, so do it now... 
        for m in self.mechanisms:
            if m not in self.ix.keys():
                raise ValueError('Mechanism %s in cell is missing from i_currents calculation', m)
 
        if showinfo:
            print('\n  [soma] find_i0  Species: %s  cell type: %s  Temp %6.1f' % (self.status['species'],
                self.status['modelType'], h.celsius))
            print('    *** found V0 = %f' % v0)
            print('    *** and cell has mechanisms: ', self.mechanisms)
        return v0



    def compute_rmrintau(self, auto_initialize=True, vrange=None):
        """
        Run the model for 2 msec after initialization - then
        compute the inverse of the sum of the conductances to get Rin at rest
        compute Cm*Rin to get tau at rest
        
        Parameters
        ----------
        auto_initialize : boolean (default: True)
            If true, forces initialization of cell in NEURON befor the computation.
            
        Returns
        -------
        A dictionary containing: Rin (Mohm), tau (ms) and Vm (mV)
        
        """
        gnames = {# R&M 03 and related:
                'nacn': 'gna', 'na': 'gna', 'jsrna': 'gna', 'nav11': 'gna', 'nacncoop': 'gna',
                'leak': 'gbar',
                'klt': 'gklt', 'kht': 'gkht',
                'ka': 'gka',
                'ihvcn': 'gh', 'hcno': 'gh', 'hcnobo': 'gh',
                # pyramidal cell specific:
                'napyr': 'gna', 'nap': 'gnap',
                'kdpyr': 'gk', 'kif': 'gkif', 'kis': 'gkis',
                'ihpyr': 'gh', 'ihpyr_adj': 'gh',
                'kcnq': 'gk',
                # cartwheel cell specific:
                'bkpkj': 'gbkpkj', 'hpkj': 'gh',
                'kpkj': 'gk', 'kpkj2': 'gk', 'kpkjslow': 'gk',
                'kpksk': 'gk', 'lkpkj': 'gbar',
                'naRsg': 'gna',
                # SGC Ih specific:
                'ihsgcApical': 'gh',  'ihsgcBasalMiddle': 'gh',
              }
        if auto_initialize:
            self.cell_initialize(vrange=vrange)
            custom_init()

        gsum = 0.
        section = self.soma
        u = self.get_mechs(section)
        for m in u:
#            gx = 'section().'+m+'.'+gnames[m]
            gm = '%s_%s' % (gnames[m], m)
            gsum += getattr(section(), gm) 
            #eval(gx)
           # print('{0:>12s} : gx '.format(m))
        # convert gsum from us/cm2 to nS using cell area
#        print ('gsum, self.somaarea: ', gsum, self.somaarea)
        gs = mho2ns(gsum, self.somaarea)
        Rin = 1e3/gs  # convert to megohms
        tau = Rin*self.totcap*1e-3  # convert to msec
        return {'Rin': Rin, 'tau': tau, 'v': self.soma(0.5).v}

    def set_soma_size_from_Cm(self, cap):
        """
        Use soma capacitance to set the cell size. Area of the open cylinder is same as a sphere of
        the same diameter.
        Compute area and save total capacitance as well
        """
        self.totcap = cap
        self.somaarea = self.totcap * 1E-6 / self.c_m  # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * ((self.somaarea / np.pi) ** 0.5)  # convert from cm to um
        self.soma.diam = lstd
        self.soma.L = lstd

    def set_soma_size_from_Diam(self, diam):
        """
        Use diameter to set the cell size. Area of the open cylinder is same as a sphere of
        the same diameter.
        Compute area and total capacitance as well
        """
        self.somaarea = 1e-8*4.*np.pi*(diam/2.)**2  # in microns^2
        self.totcap = self.c_m * self.somaarea * 1e6
    #    lstd = diam # 1E4 * ((self.somaarea / np.pi) ** 0.5)  # convert from cm to um
        self.soma.diam = diam
        self.soma.L = diam
        
    def set_soma_size_from_Section(self, soma):
        self.soma.diam = soma.diam
        self.soma.L = soma.L
        self.somaarea = 1e-8*np.pi*soma.diam*soma.L
        self.totcap = self.c_m * self.somaarea * 1e6

    def print_soma_info(self):
        print('-'*40)
        print('Soma Parameters: ')
        print('   Area: ', self.somaarea)
        print('   Cap:  ', self.totcap)
        print('   L:    ', self.soma.L)
        print('   diam: ', self.soma.diam)
        print('   cm:   ', self.c_m)
        print('-'*40)
        
    def distances(self, section=None):
        self.distanceMap = {}
        if section is None:
            self.hr.h('access %s' % self.soma.name()) # reference point
        else:
            self.hr.h('access %s' % section.name())
        d = self.hr.h.distance()
        for sec in self.all_sections:
            s = self.all_sections[sec]
            if len(s) > 0:
                for u in s:
                    self.hr.h('access %s' % u.name())
                    self.distanceMap[u.name()] = self.hr.h.distance(0.5) - d # should be distance from first point

    def computeAreas(self):
        self.areaMap = {}
        for sec in self.all_sections:  # keys for names of section types
            s = self.all_sections[sec]  # get all the sections of that type
            sectype = self.get_section_type(s)
            if len(s) > 0:
                self.areaMap[sec] = {}
                for u in s:
                    self.areaMap[sec][u] = np.pi*u.diam*u.L
            else:
                pass
                # print('    No section of type %s in cell' % sec)

    def add_axon(self, c_m=1.0, R_a=150, axonsf=1.0, nodes=5, debug=False, dia=None, len=None, seg=None):
        """
        Add an axon to the soma with an initial segment (tapered), and multiple nodes of Ranvier
        The size of the axon is determined by self.axonsf, which in turn is set by the species
        The somaarea is used to scale the density of ion channels in the initial segment
        """
        nnodes = range(nodes)
        axnode = []
        internode = []
        Section = h.Section
        initsegment = Section(cell=self.soma)
        initsegment.connect(self.soma)
        for i in nnodes:
            axnode.append(Section(cell=self.soma))
            internode.append(Section(cell=self.soma))
        axnode[0].connect(initsegment)
        for i in nnodes:
            internode[i].connect(axnode[i])
            if i < nnodes[-1]:
                axnode[i + 1].connect(internode[i])

                # create an initial segment
        ninitseg = 21
        initsegment.nseg = ninitseg
        initsegment.diam = 4.0 * axonsf
        initsegment.L = 36.0 * axonsf
        initsegment.cm = c_m # c_m
        initsegment.Ra = R_a # R_a
        initsegment.insert('nacn')  # uses a standard Rothman sodium channel
        initsegment.insert('kht')
        initsegment.insert('klt')
        initsegment.insert('ihvcn')
        initsegment.insert('leak')
        gnamax = nstomho(6000.0, self.somaarea)
        gnamin = 0.0 * gnamax

        gnastep = (gnamax - gnamin) / ninitseg  # taper sodium channel density
        for ip, inseg in enumerate(initsegment):
            gna = gnamin + ip * gnastep
            if debug:
                print('Initial segment %d: gnabar = %9.6f' % (ip, gna))
            inseg.nacn.gbar = gna
            inseg.klt.gbar = 0.2 * nstomho(200.0, self.somaarea)
            inseg.kht.gbar = nstomho(150.0, self.somaarea)
            inseg.ihvcn.gbar = 0.0 * nstomho(20.0, self.somaarea)
            inseg.leak.gbar = nstomho(2.0, self.somaarea)
            inseg.ena = self.e_na
            inseg.ek = self.e_k
            inseg.leak.erev = self.e_leak

        for i in nnodes:
            axnode[i] = self.loadaxnodes(axnode[i], self.somaarea, eleak=self.e_leak)
            internode[i] = self.loadinternodes(internode[i], self.somaarea, eleak=self.e_leak)

        if debug:
            print("<< {:s} Axon Added >>".format(self.__class__.__name__))
            h.topology()
        self.add_section(initsegment, 'initialsegment')
        self.add_section(axnode, 'axonnode')
        self.add_section(internode, 'internode')

    @staticmethod
    def loadaxnodes(axnode, somaarea, nodeLength=2.5, nodeDiameter=2.0, eleak=-65):
        v_potassium = -80  # potassium reversal potential
        v_sodium = 50  # sodium reversal potential
        Ra = 150
        cm = 1.0
        axnode.nseg = 1
        axnode.L = nodeLength
        axnode.diam = nodeDiameter
        axnode.Ra = Ra
        axnode.cm = cm
        axnode.insert('nacn')
        axnode.insert('kht')
        axnode.insert('klt')
        axnode.insert('leak')
        axnode.insert('ihvcn')
        for ax in axnode:
            ax.nacn.gbar = nstomho(1000.0, somaarea)
            ax.kht.gbar = nstomho(150.0, somaarea)
            ax.klt.gbar = nstomho(200.0, somaarea)
            ax.ihvcn.gbar = 0
            ax.leak.gbar = nstomho(2.0, somaarea)
            ax.ena = v_sodium
            ax.ek = v_potassium
            ax.leak.erev = eleak
        return axnode

    @staticmethod
    def loadinternodes(internode, somaarea, internodeLength=1000, internodeDiameter=10, eleak=-65):
        v_potassium = -80  # potassium reversal potential
        v_sodium = 50  # sodium reversal potential
        Ra = 150
        cm = 0.002

        internode.nseg = 20
        internode.L = internodeLength
        internode.diam = internodeDiameter
        internode.Ra = Ra
        internode.cm = cm
        internode.insert('nacn')
        internode.insert('kht')
        internode.insert('leak')
        for inno in internode:
            inno.leak.gbar = nstomho(0.002, somaarea)
            inno.nacn.gbar = 0 * nstomho(500.0, somaarea)
            inno.kht.gbar = 0 * nstomho(150.0, somaarea)
            inno.ek = v_potassium
            inno.ena = v_sodium
            inno.leak.erev = eleak
        return internode
