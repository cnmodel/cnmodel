Things tested and to do:

ccstim: called by test_stim, but does ot appear to exist. 
stim does not provide all the protocols tested in test_stim

test_circuit.py:
      File "test_circuit.py", line 42, in <module>
        bushy.resolve_inputs(depth=2)
      File "/Users/pbmanis/Desktop/Python/cnmodel/cnmodel/populations/population.py", line 133, in resolve_inputs
        pre_cells = self.connect_pop_to_cell(pop, i)
      File "/Users/pbmanis/Desktop/Python/cnmodel/cnmodel/populations/population.py", line 166, in connect_pop_to_cell
        pre_cell.connect(cell)
      File "/Users/pbmanis/Desktop/Python/cnmodel/cnmodel/cells/cell.py", line 217, in connect
        synapse = synapses.Synapse(self, pre_opts, post_cell, post_opts, **kwds)
      File "/Users/pbmanis/Desktop/Python/cnmodel/cnmodel/synapses/synapse.py", line 13, in __init__
        self.terminal = pre_cell.make_terminal(post_cell, **pre_opts)
      File "/Users/pbmanis/Desktop/Python/cnmodel/cnmodel/cells/cell.py", line 246, in make_terminal
        post_cell.__class__.__name__))
    NotImplementedError: Cannot make Terminal connecting TStellateRothman => BushyRothman

test_mechanisms:
    some mechanisms do not make sense to be tested in this routine; need to exclude (done)
    
    
test_physiology: ok

test_populations : not sure what this is supposed to do.

test_sgc_input.py : works.

test_sgc_input_phaselocking: no bushy cell spikes. Is inputs strength scaled correctly?

test_simple_synapses works. 

test_sound_stim works

test_synapses (more detailed test) works

toy_model : works.





