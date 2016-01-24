cnmodel
=======

cnmodel is a Python-based interface to NEURON models of cochlear nucleus neurons, synapses, network connectivity, and the Zilaney et al model of auditory nerve fiber spike trains to acoustic stimuli. The overall goal is to provide a platform for modeling networks of cochlear nucleus neurons with different levels of biological realism in the context of an accurate simulation of auditory nerve input.

At the lowest level are NEURON-based implementations of ion channels and synapses. Ion channel models for potassium channels are derived largely from the measurements and models of Rothman and Manis (2003abc), and Kanold and Manis (1999, 2001); other channels are derived or modified from the literature. Cell models are in turn based on the insertion of ion channels in densities based on measurements in guinea pig and mouse. The "point" somatic cell models, which form the base set of models in cnmodel, replicate the data reported in the original papers. 

The postsynaptic receptor/conductance models are based on kinetic models of glutamate (Raman and Trussell, 1992) and glycinergic receptors, as adjusted to match measurements of synaptic conductances from the mouse cochlear nucleus collected in Xie and Manis, 2013. The glutamate receptor models include desensitization and the effects of internal polyamine receptor block, based on the kinetic scheme of Woodhull (1982).

The presynaptic release model includes a multisite, probabilistic synapse that includes time-dependent changes in release probability based on the Dittman, Kreiger and Regehr (2000) kinetic scheme. Although detailed, this model is computationally expensive and not suitable for large scale network simulations. 

Network connectivity may be defined programatically, or based on a table of connetivity patterns. A table with estimates derived from the literature is included in the source. 

Included in this package is a set of test suites for different components. An overriding set of unit tests is available to confirm performance of the base models against a set of reference runs, several of which are in turn directly traceable to the original manuscripts. The test suites are useful in verifying performance of the model after modifications of the code or changes in the installation (upgrades of Python or Matlab, for example). 

Installation requirements
-------------------------
This package depends on the following:

   1. Python 2.7.10 with at least numpy, scipy, Qt4.8 and pyqtgraph. An Anaconda install with the appropriate scientific packages works well.
   2. A Python-linked version of NEURON (www.neuron.yale.edu). The code has been tested with NEURON 7.3 and 7.4.
   3. A C compiler (gcc). Needed for compilation of C files for matlab, and compilation of the mechanisms for NEURON.
   4. Matlab (2011 or later). Matlab is needed for the auditory nerve model to run. The C compiler is needed to compile the C files in the auditory nerve model. 
   5. neuronvis (available at https://github.com/campagnola/neuronvis). This provides a library that can read and visualize hoc files.

It is easiest to set up a model using a separate directory, with links (ln -s) to cnmodel and neuronvis, rather than installing those as libraries in the Python environment.

Testing
-------
Before testing, enter the cnmodel directory, and compile the NEURON mod files:

    nrnivmodl cnmodel/mechanisms

This will create a directory ("x86_64" or "special") in the top cnmodel directory. At that point

    python toy_model.py
     
should generate a plot of 9 sets of traces showing responses of individual neuron models to depolarizing and hyperpolarizing current steps.
The test suite can be run as:

    python test.py

This will test each of the models against reference data, the synapse mechanisms, a number of internal routines, and the auditory nerve model. The tests should pass for each component. Failures may indicate incorrect installation or functioning of individual components.



