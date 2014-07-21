"""
Wrapper classes to make working with NEURON easier.

Author: Andrew P. Davison, UNIC, CNRS
"""

__version__ = "0.3.0"

from neuron import nrn, h, hclass

h.load_file('stdrun.hoc')

PROXIMAL = 0
DISTAL = 1

class Mechanism(object):
    """
    Examples:
    >>> leak = Mechanism('pas', {'e': -65, 'g': 0.0002})
    >>> hh = Mechanism('hh')
    added set_parameters to allow post-instantiation parameter modification
    """
    def __init__(self, name, **parameters):
        self.name = name
        self.parameters = parameters

    def set_parameters(self, parameters):
        self.parameters = parameters

    def insert_into(self, section):
        section.insert(self.name)
        for name, value in self.parameters.items():
            for segment in section:
                mech = getattr(segment, self.name)
                setattr(mech, name, value)


class Section(nrn.Section):
    """
    Examples:
    >>> soma = Section(L=30, diam=30, mechanisms=[hh, leak])
    >>> apical = Section(L=600, diam=2, nseg=5, mechanisms=[leak],
    ...                  parent=soma, connection_point=DISTAL)
    """
    
    def __init__(self, L, diam, nseg=1, Ra=100, cm=1, mechanisms=[], parent=None, connection_point=DISTAL):
        nrn.Section.__init__(self)
        # set geometry
        self.L = L
        self.diam = diam
        self.nseg = nseg
        # set cable properties
        self.Ra = Ra
        self.cm = cm
        # connect to parent section
        if parent:
            self.connect(parent, connection_point, PROXIMAL)
        # add ion channels
        for mechanism in mechanisms:
            mechanism.insert_into(self)

    def add_synapses(self, label, type, locations=[0.5], **parameters):
        if hasattr(self, label):
            raise Exception("Can't overwrite synapse labels (to keep things simple)")
        synapse_group = []
        for location in locations:        
            synapse = getattr(h, type)(location, sec=self)
            for name, value in parameters.items():
                setattr(synapse, name, value)
            synapse_group.append(synapse)
        if len(synapse_group) == 1:
            synapse_group = synapse_group[0]
        setattr(self, label, synapse_group)
    add_synapse = add_synapses  # for backwards compatibility

    def plot(self, variable, location=0.5, tmin=0, tmax=5, xmin=-80, xmax=40):
        import neuron.gui
        self.graph = h.Graph()
        h.graphList[0].append(self.graph)
        self.graph.size(tmin, tmax, xmin, xmax)
        self.graph.addvar('%s(%g)' % (variable, location), sec=self)

    def record_spikes(self, threshold=-30):
        self.spiketimes = h.Vector()
        self.spikecount = h.APCount(0.5, sec=self)
        self.spikecount.thresh = threshold
        self.spikecount.record(self.spiketimes)


def alias(attribute_path):
    """
    Returns a new property, mapping an attribute nested in an object hierarchy
    to a simpler name

    For example, suppose that an object of class A has an attribute b which
    itself has an attribute c which itself has an attribute d. Then placing
      e = alias('b.c.d')
    in the class definition of A makes A.e an alias for A.b.c.d
    """
    parts = attribute_path.split('.')
    attr_name = parts[-1]
    attr_path = parts[:-1]
    def set(self, value):
        obj = reduce(getattr, [self] + attr_path)
        setattr(obj, attr_name, value)
    def get(self):
        obj = reduce(getattr, [self] + attr_path)
        return getattr(obj, attr_name)
    return property(fset=set, fget=get)


def uniform_property(section_list, attribute_path):
    """
    Define a property that will have a uniform value across a list of sections.
    
    For example, suppose we define a neuron model as a class A, which contains
    three compartments: soma, dendrite and axon. Then placing
    
        gnabar = uniform_property(["soma", "axon"], "hh.gnabar")
    
    in the class definition of A means that setting a.gnabar (where a is an
    instance of A) will set the value of hh.gnabar in both the soma and axon, i.e.

        a.gnabar = 0.01
        
    is equivalent to:
    
        for sec in [a.soma, a.axon]:
            for seg in sec:
                seg.hh.gnabar = 0.01

    """
    parts = attribute_path.split('.')
    attr_name = parts[-1]
    attr_path = parts[:-1]
    def set(self, value):
        for sec_name in section_list:
            sec = getattr(self, sec_name)
            for seg in sec:
                obj = reduce(getattr, [seg] + attr_path)
                setattr(obj, attr_name, value)
    def get(self):
        sec = getattr(self, section_list[0])
        obj = reduce(getattr, [sec(0.5)] + attr_path)
        return getattr(obj, attr_name)
    return property(fset=set, fget=get)



if __name__ == "__main__":
    
    class SimpleNeuron(object):
    
        def __init__(self):
            # define ion channel parameters
            leak = Mechanism('pas', e=-65, g=0.0002)
            hh = Mechanism('hh')
            # create cable sections
            self.soma = Section(L=30, diam=30, mechanisms=[hh])
            self.apical = Section(L=600, diam=2, nseg=5, mechanisms=[leak], parent=self.soma,
                                  connection_point=DISTAL)
            self.basilar = Section(L=600, diam=2, nseg=5, mechanisms=[leak], parent=self.soma,
                                   connection_point=0.5)
            self.axon = Section(L=1000, diam=1, nseg=37, mechanisms=[hh],
                                connection_point=0)
            # synaptic input
            self.soma.add_synapses('syn', 'AlphaSynapse', onset=0.5, gmax=0.05, e=0)
    
        gnabar = uniform_property(["soma", "axon"], "hh.gnabar")
        gkbar = uniform_property(["soma", "axon"], "hh.gkbar")
    
    neuron = SimpleNeuron()
    neuron.soma.plot('v')
    neuron.apical.plot('v')
    
    print neuron.gnabar
    neuron.gnabar = 0.15
    assert neuron.soma(0.5).hh.gnabar == 0.15


    h.dt = 0.025
    v_init = -65
    tstop = 5
    h.finitialize(v_init)
    h.run()
