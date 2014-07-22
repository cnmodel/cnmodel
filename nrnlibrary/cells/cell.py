from neuron import h
import scipy.optimize
import numpy as np
from ..util import nstomho, mho2ns

class Cell(object):
    """
    Base class for all cell types.
    """
    sec_lookup = {}  # create a lookup table to map sections to their parent cell
    
    def __init__(self):
        # dictionary of all sections associated with this cell
        self.all_sections = {}
        # the following section types (parts) are known to us:
        for k in ['soma', 'maindend', 'secdend', 'internode', 'initialsegment', 'axonnode']:
            self.all_sections[k] = []  # initialize to an empty list
        
        # each cell has the following parameters:
        self.totcap = None  # total membrane capacitance (somatic)
        self.somaarea = None  # total soma area
        self.initsegment = None  # hold initial segment sections
        self.axnode = None  # hold nodes of ranvier sections
        self.internode = None  # hold internode sections
        self.maindend = None  # hold main dendrite sections
        self.secdend = None  # hold secondary dendrite sections
        self.axonsf = None  # axon diameter scale factor
        # define defaults for these parameters (RM03 model defaults)
        self.e_k = -70  # potassium reversal potential, mV
        self.e_na = 55
        self.e_h = -43
        self.c_m = 0.9  # specific membrane capacitance,  uf/cm^2
        self.R_a = 150  # axial resistivity of cytoplasm/axoplasm, ohm.cm
        self.e_leak = -65

        # Recommended current (min, max, step) for testing this cell
        self.i_test_range=(-0.5, 0.5, 0.05)
        
        # Recommended threshold for detecting spikes from this cell
        self.spike_threshold = -40
        
        # Resting potential for this cell, determined by calling
        # self.find_i0()
        self.vm0 = None

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
    
    @property
    def soma(self):
        """
        First (or only) section in the "soma" section group.
        """
        return self.all_sections['soma'][0]

    def print_status(self):
        print("\nCell model: %s" % self.__class__.__name__)
        print(self.__doc__)
        print '    Model Status:'
        print '-'*24
        for s in self.status.keys():
            print('{0:>12s} : {1:<12s}'.format(s, repr(self.status[s])))
        print '-'*32

    def cell_initialize(self, showinfo=False):
        """
        Initialize this cell to it's "rmp" under current conditions
        All sections in the cell are set to the same value
        """
        if self.vm0 is None:
            self.vm0 = self.find_i0(showinfo=showinfo)
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = self.vm0

    def get_mechs(self, section):
        """
        return a list of the mechanisms that are present in a section
        a mechanism is required to have a gbar variable.
        """
        u=dir(section())
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
        print '\n    Installed mechanisms:'
        self.get_mechs(section)
        for m in self.mechs:
            try:
                gx=eval('section().'+m+'.gbar')
                print gx
                print('{0:>12s} : {1:<7.3g} mho/cm2  {2:<7.3g} nS '.format(m, repr(gx), repr(mho2ns(gx))))
            except:
                print('{0:>12s} : <no gbar> '.format(m))
        print '-'*32


    def i_currents(self, V):
        """
        For the steady-state case, return the total current at voltage V
        Used to find the zero current point
        vrange brackets the interval

        """
        self.ix = {}

        # klt
        if 'klt' in self.mechanisms:
            zss = 0.5   # <0,1>   : steady state inactivation of glt
            q10g = 2.0
            qg = q10g**((h.celsius-22.)/10.)
            winf = (1.0 / (1.0 + np.exp(-(V + 48.) / 6 )))**0.25
            zinf = zss + ((1.0-zss) / (1.0 + np.exp((V + 71.) / 10 )))
            gklt = qg*self.soma().klt.gbar*(winf**4.)*zinf
            self.ix['klt'] = gklt*(V - self.soma().ek)

        # ihvcn
        if 'ihvcn' in self.mechanisms:
            q10g = 2.0
            qg = q10g**((h.celsius-22.)/10.)
            rinf = 1. / (1+np.exp((V + 76.) / 7.))
            gh = qg*self.soma().ihvcn.gbar*rinf
            self.ix['ihvcn'] = gh*(V - self.soma().ihvcn.eh)

        # kht
        if 'kht' in self.mechanisms:
            q10g = 2.0
            qg = q10g**((h.celsius-22)/10)
            nf = 0.85  # proportion of n vs p kinetics
            ninf =   (1 + np.exp(-(V + 15.) / 5.))**(-0.5)
            pinf =  1 / (1 + np.exp(-(V + 23.) / 6.))
            gkht = qg*self.soma().kht.gbar*(nf*(ninf**2) + (1.-nf)*pinf)
            self.ix['kht'] = gkht*(V - self.soma().ek)

        # ka
        if 'ka' in self.mechanisms:
            q10g = 2.0
            qg = q10g**((h.celsius-22)/10)
            ainf = (1. / (1 + np.exp(-1*(V + 31.) / 6 )))**0.25
            binf = 1 / (1 + np.exp((V + 66.) / 7))**0.5
            cinf = 1 / (1 + np.exp((V + 66.) / 7))**0.5
            gka = self.soma().ka.gbar*(ainf**4)*binf*cinf
            self.ix['ka'] = gka*(V - self.soma().ek)

        # na
        if 'na' in self.mechanisms:
            q10g = 2.0
            qg = q10g**((h.celsius-22.)/10.)
            minf = 1 / (1+np.exp(-(V + 38.) / 7.))
            hinf = 1 / (1+np.exp((V + 65.) / 6.))
            gna = qg*self.soma().na.gbar*(minf**3.0)*hinf
            self.ix['na'] = gna*(V - self.soma.ena)

        # leak
        if 'leak' in self.mechanisms:
            self.ix['leak'] = self.soma().leak.gbar*(V - self.soma().leak.erev)
    #    print ix
        return np.sum([self.ix[i] for i in self.ix])

    def find_i0(self, vrange=[-70., -55.], showinfo=False):
        """
        find the root of the system of equations in vrange.
        Finds RMP fairly accurately as zero current level for current conductances.
        """
        v0 = scipy.optimize.brentq(self.i_currents, vrange[0], vrange[1])
        if showinfo:
            print '\n  find_i0  Species: %s  cell type: %s' % (self.status['species'], self.status['type'])
            print '    *** found V0 = %f' % v0
            print '    *** using conductances: ', self.ix.keys()
            print '    *** and cell has mechanisms: ', self.mechanisms
        return v0

    def measure_rintau(self, auto_initialize = True):
        """
        Run the model for 2 msec after initialization - then
        compute the inverse of the sum of the conductances to get Rin at rest
        compute Cm*Rin to get tau at rest
        :param none:
        :return Rin (Mohm), tau (ms) and Vm (mV):
        """
        if auto_initialize:
            self.cell_initialize()
        gnames = {'nacn': 'gna', 'na': 'gna',
                  'leak': 'gbar',
                  'klt': 'gklt', 'kht': 'gkht',
                  'ka': 'gka',
                  'ihvcn': 'gh', 'hcno': 'thegna'
                  }
        gsum = 0.
        section = self.soma
        u = self.get_mechs(section)
        for m in u:
            gx = 'section().'+m+'.'+gnames[m]
            gsum += eval(gx)
           # print('{0:>12s} : gx '.format(m))
        # convert gsum from us/cm2 to nS using cell area
        gs = mho2ns(gsum, self.somaarea)
        Rin = 1e3/gs  # convert to megohms
        tau = Rin*self.totcap*1e-3  # convert to msec
        return Rin, tau, self.soma(0.5).v

    def set_soma_size_from_Cm(self, cap):
        self.totcap = cap
        self.somaarea = self.totcap * 1E-6 / self.c_m  # pf -> uF, cm = 1uf/cm^2 nominal
        lstd = 1E4 * ((self.somaarea / np.pi) ** 0.5)  # convert from cm to um
        self.soma.diam = lstd
        self.soma.L = lstd

    def add_axon(self, c_m=1.0, R_a=150, axonsf=1.0, nodes=5, debug=False):
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
        initsegment.cm = c_m
        initsegment.Ra = R_a
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
                print 'Initial segment %d: gnabar = %9.6f' % (ip, gna)
            inseg.nacn.gbar = gna
            inseg.klt.gbar = 0.2 * nstomho(200.0, self.somaarea)
            inseg.kht.gbar = nstomho(150.0, self.somaarea)
            inseg.ihvcn.gbar = 0.0 * nstomho(20.0, self.somaarea)
            inseg.leak.gbar = nstomho(2.0, self.somaarea)
            inseg.ena = self.e_na
            inseg.ek = self.e_k

        for i in nnodes:
            axnode[i] = self.loadaxnodes(axnode[i], self.somaarea)
            internode[i] = self.loadinternodes(internode[i], self.somaarea)

        if debug:
            print("<< {:s} Axon Added >>".format(self.__class__.__name__))
            h.topology()
        self.add_section(initsegment, 'initsegment')
        self.add_section(axnode, 'axonnode')
        self.add_section(internode, 'internode')

    @staticmethod
    def loadaxnodes(axnode, somaarea, nodeLength=2.5, nodeDiameter=2.0):
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
        return axnode

    @staticmethod
    def loadinternodes(internode, somaarea, internodeLength=1000, internodeDiameter=10):
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
            inno.leak.e = -80
        return internode
