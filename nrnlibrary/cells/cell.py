from neuron import h
from ..pynrnutilities import nstomho

class Cell(object):
    """
    Base class for all cell types.
    """
    def __init__(self):
        pass

    def print_status(self):
        print("\nCell model: %s" % self.__class__.__name__)
        print(self.__doc__)
        print '    Model Status:'
        print '-'*24
        for s in self.status.keys():
            print('{0:>12s} : {1:<12s}'.format(s, repr(self.status[s])))
        print '-'*32

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
                print('{0:>12s} : {1:<7.3g} '.format(m, eval('section().'+m+'.gbar')))
            except:
                print('{0:>12s} : <no gbar> '.format(m))
        print '-'*32

    def add_axon(self, soma, somaarea, c_m=1.0, R_a=150, axonsf=1.0, nodes=5, debug=False):
        """
        Add an axon to the soma with an initial segment (tapered), and multiple nodes of Ranvier
        The size of the axon is determined by self.axonsf, which in turn is set by the species
        The somaarea is used to scale the density of ion channels in the initial segment
        """
        nnodes = range(nodes)
        axnode = []
        internode = []
        Section = h.Section
        initsegment = Section(cell=soma)
        initsegment.connect(soma)
        for i in nnodes:
            axnode.append(Section(cell=soma))
            internode.append(Section(cell=soma))
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
        gnamax = nstomho(6000.0, somaarea)
        gnamin = 0.0 * gnamax

        gnastep = (gnamax - gnamin) / ninitseg  # taper sodium channel density
        for ip, inseg in enumerate(initsegment):
            gna = gnamin + ip * gnastep
            if debug:
                print 'Initial segment %d: gnabar = %9.6f' % (ip, gna)
            inseg.nacn.gbar = gna
            inseg.klt.gbar = 0.2 * nstomho(200.0, somaarea)
            inseg.kht.gbar = nstomho(150.0, somaarea)
            inseg.ihvcn.gbar = 0.0 * nstomho(20.0, somaarea)
            inseg.leak.gbar = nstomho(2.0, somaarea)
            inseg.ena = self.e_na
            inseg.ek = self.e_k

        for i in nnodes:
            axnode[i] = self.loadaxnodes(axnode[i], self.somaarea)
            internode[i] = self.loadinternodes(internode[i], self.somaarea)

        if debug:
            print("<< {:s} Axon Added >>".format(self.__class__.__name__))
            h.topology()
        return(initsegment, axnode, internode)

    @staticmethod
    def loadaxnodes(axnode, somaarea, scalefactor, nodeLength=2.5, nodeDiameter=2.0):
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
    def loadinternodes(internode, somaarea, scalefactor, internodeLength=1000, internodeDiameter=10):
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
