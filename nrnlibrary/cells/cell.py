from neuron import h
from ..util import nstomho, mho2ns

class Cell(object):
    """
    Base class for all cell types.
    """
    sec_lookup = {}  # create a lookup table to map sections to their parent cell
    def add_section(self, sec):
        Cell.sec_lookup[sec.name()] = self

    def __init__(self):
        self.all_sections = {}  # dictionary of all sections associated with this cell
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

        # Recommended current (min, max, step) for testing this cell
        self.i_test_range=(-0.5, 0.5, 0.05)
        
        # Recommended threshold for detecting spikes from this cell
        self.spike_threshold = -40


    def print_status(self):
        print("\nCell model: %s" % self.__class__.__name__)
        print(self.__doc__)
        print '    Model Status:'
        print '-'*24
        for s in self.status.keys():
            print('{0:>12s} : {1:<12s}'.format(s, repr(self.status[s])))
        print '-'*32

    def initialize(self):
        """
        Initialize this cell to it's specified RMP.
        All sections in the cell are set to the same value
        """
        for part in self.all_sections.keys():
            for sec in self.all_sections[part]:
                sec.v = self.vm0
        h.finitialize()
        h.fcurrent()


    def print_mechs(self, section):
        """
        print the mechanisms that are inserted into the specified section,
        and their densities (in uS/cm^2)
        """
        print '\n    Installed mechanisms:'
        u=dir(section())
        for m in u:
            if m[0:2] == '__':
                continue
            if m in ['cm', 'diam', 'k_ion', 'na_ion', 'next', 'point_processes', 'sec', 'v', 'x']:
                continue  # skip non-mechanisms known to us
            try:
                gx=eval('section().'+m+'.gbar')
                print gx
                print('{0:>12s} : {1:<7.3g} mho/cm2  {2:<7.3g} nS '.format(m, repr(gx), repr(mho2ns(gx))))
            except:
                print('{0:>12s} : <no gbar> '.format(m))
        print '-'*32

    def measure_rintau(self):
        """
        Run the model for 2 msec after initialization - then
        compute the inverse of the sum of the conductances to get Rin at rest
        compute Cm*Rin to get tau at rest
        :param none:
        :return Rin (Mohm), tau (ms) and Vm (mV):
        """
        u=dir(self.soma())
        h.tstop = 2.0
        h.init()
        h.finitialize(self.vm0)
        h.run()
        gnames = {'nacn': 'gna', 'leak': 'gbar',
                  'klt': 'gklt', 'kht': 'gkht',
                  'ka': 'gka','ihvcn': 'gh',
                  }
        gsum = 0.
        for m in u:
            if m[0:2] == '__':
                continue
            if m in ['cm', 'diam', 'k_ion', 'na_ion', 'next', 'point_processes', 'sec', 'v', 'x']:
                continue  # skip non-mechanisms known to us
            gx = 'section().'+m+'.'+gnames[m]
            try:
                gsum += eval(gx)
            except:
                pass
                #print('{0:>12s} : <no g> '.format(m))
        # convert gsum from us/cm2 to nS using cell area
        gs = mho2ns(gsum, self.somaarea)
        Rin = 1e3/gs  # convert to megohms
        tau = Rin*self.totcap*1e-3  # convert to msec
       # print('Vm0: {0:8.1f} mV     Vmeas: {1:8.2f} mV'.format(self.vm0, self.soma(0.5).v))
       # print('Rin: {0:8.1f} Mohm   tau: {1:8.2f} ms  cap: {2:8.1f}'.format(Rin, tau, self.totcap))
        return Rin, tau, self.som(0.5).v

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
        self.all_sections['initsegment'].extend(initsegment)
        self.all_sections['axonnode'].extend(axnode)
        self.all_sections['internode'].extend(internode)

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
