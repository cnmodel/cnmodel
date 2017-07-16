"""
Test construction of a complete circuit build from cell populations.

This script:

1. Creates populations of sgc, bushy, and stellate cells.
2. Connects the populations together.
3. Instantiates 10 bushy and 10 t-stellate cells near 16kHz.
4. Resolves all required synaptic dependencies.

No simulation is executed; this is meant to be run interactively to allow 
introspection of the circuit that was generated.
"""

from cnmodel import populations

def testcircuit():
    # Create cell populations.
    # This creates a complete set of _virtual_ cells for each population. No 
    # cells are instantiated at this point.
    sgc = populations.SGC()
    bushy = populations.Bushy()
    dstellate = populations.DStellate()
    tstellate = populations.TStellate()

    # Connect populations. 
    # This only records the connections between populations; no synapses are 
    # created at this stage.
    sgc.connect(bushy, tstellate, dstellate)
    dstellate.connect(tstellate, bushy)
    #tstellate.connect(bushy)  # this will fail - we don't know about this connection yet.

    # Select cells to record from.
    # At this time, we actually instantiate the selected cells.
    # select 10 bushy cells closest to 16kHz
    bushy_cell_ids = bushy.select(10, cf=16e3, create=True)  
    # select 10 stellate cells closest to 16kHz
    tstel_cell_ids = tstellate.select(10, cf=16e3, create=True)  

    # Now create the supporting circuitry needed to drive the cells we selected.
    # At this time, cells are created in all populations and automatically 
    # connected with synapses.
    bushy.resolve_inputs(depth=2)
    tstellate.resolve_inputs(depth=2)
    # Note that using depth=2 indicates the level of recursion to use when 
    # resolving inputs. For example, resolving inputs for the bushy cell population
    # (level 1) creates presynaptic cells in the dstellate population, and resolving
    # inputs for the dstellate population (level 2) creates presynaptic cells in the
    # sgc population. 


    # TODO: 
    #   - specify which parameters to record (Vm, spike times, per-synapse currents, etc)
    #   - run simulation and display / analyze results
    #   - add examples of modifying circuitry to search parameter spaces, test
    #     hypotheses, etc.
    
if __name__ == '__main__':
    testcircuit()