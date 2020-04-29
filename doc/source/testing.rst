Testing
=======

Make sure you are in the cnmodel directory, and that you have selected the right environment in Anaconda (in 
my case, this is usually an environment called py3mpl3 - python 3 with matplotlib 3).

After the code is installed, enter the cnmodel directory and compile the NEURON mod (you might have already done this if you are following the instructions above)::

    $ nrnivmodl cnmodel/mechanisms

This will create a directory ("x86_64" or "special") in the top cnmodel directory with the compiled mechanisms.

Then::

    $ python examples/toy_model.py
     
should generate a plot with several sets of traces showing responses of individual neuron models to depolarizing and hyperpolarizing current steps.

The test suite should then be run as::

    $ python test.py

This will test each of the models against reference data, the synapse mechanisms, a number of internal routines, and the auditory nerve model. The tests should pass for each component. Failures may indicate incorrect installation or incorrect function within individual components. These should be corrected before proceeding with simulations.

Individual test suite components can be run directly using pytest, for example::

    $  pytest cnmodel/cells/tests/test_cells.py -k "test_bushy_mouse"