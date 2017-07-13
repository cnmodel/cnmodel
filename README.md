About CNModel
=============

CNModel is a Python-based interface to NEURON models of cochlear nucleus neurons, synapses, network connectivity. To drive the model with sound stimuli, a bridge to the MATLAB model of Zilaney et al (2010, 2014)  of auditory nerve fiber spike trains to acoustic stimuli is provided, or a Python-based version of the model from Rudnicki and Hemmert can be used. The overall goal is to provide a platform for modeling networks of cochlear nucleus neurons with different levels of biological realism in the context of an accurate simulation of auditory nerve input.

At the lowest level are NEURON-based implementations of ion channels and synapses. Ion channel models for potassium channels are derived largely from the measurements and models of Rothman and Manis (2003abc), and Kanold and Manis (1999, 2001); other channels are derived or modified from the literature. Cell models are in turn based on the insertion of ion channels in densities based on measurements in guinea pig and mouse. The "point" somatic cell models, which form the base set of models in CNModel, replicate the data reported in the original papers. 

The postsynaptic receptor/conductance models are based on kinetic models of glutamate (Raman and Trussell, 1992) and glycinergic receptors, as adjusted to match measurements of synaptic conductances from the mouse cochlear nucleus collected in Xie and Manis, 2013. The glutamate receptor models include desensitization and the effects of internal polyamine receptor block, based on the kinetic scheme of Woodhull (1982).

The presynaptic release model includes a multisite, probabilistic synapse that includes time-dependent changes in release probability based on the Dittman, Kreitzer and Regehr (J Neurosci. 2000 Feb 15;20(4):1374-85) kinetic scheme. Although detailed, this model is computationally expensive and likely not suitable for large scale network simulations. Other simpler models of synapses are also included.

Network connectivity may be defined programmatically, or based on a table of connectivity patterns. A table with estimates derived from the literature is included in the source. 

Included in this package is a set of test suites for different components. An overriding set of unit tests is available to confirm performance of the base models against a set of reference runs, several of which are in turn directly traceable to the original manuscripts. The test suites are useful in verifying performance of the model after modifications of the code or changes in the installation (upgrades of Python or Matlab, for example). 

Installation requirements
-------------------------
This package depends on the following:

   1. Python 2.7.10 with numpy (1.11), scipy (0.19.0), lmfit (0.9.6), Qt4.8.7, pyqt4 (4.11.4) and pyqtgraph (0.9.10). An Anaconda install with the appropriate scientific packages works well. lmfit is best obtained via pip to install the latest versions. This package is not yet compatible with Python 3.x.
   2. neuronvis (available at https://github.com/campagnola/neuronvis or https://github.com/pbmanis/neuronvis). This provides a library that can read and visualize hoc files.
   3. A Python-linked version of NEURON (www.neuron.yale.edu). The code has been tested with NEURON 7.3 and 7.4.
   4. A C compiler (gcc). Needed for compilation of C files for matlab, and compilation of the mechanisms for NEURON.
   5. The Zilany et al (JASA 2014) auditory periphery model. This can be provided one of two ways:
      1. The original MATLAB-based Zilany model; requires MATLAB 2011 or later. A C compiler will also be needed to build this model.
      2. The Python-based cochlea model (https://github.com/mrkrd/cochlea)   5. The Zilany et al (JASA 2014) cochlea model (available at https://www.urmc.rochester.edu/MediaLibraries/URMCMedia/labs/carney-lab/codes/Zilany-2014-Code-and-paper.zip).


Testing
-------
Before testing, enter the CNModel directory, and compile the NEURON mod files:

    nrnivmodl CNModel/mechanisms

This will create a directory ("x86_64" or "special") in the top CNModel directory. At that point

    python toy_model.py
     
should generate a plot with several sets of traces showing responses of individual neuron models to depolarizing and hyperpolarizing current steps.


The test suite can be run as:

    python test.py

This will test each of the models against reference data, the synapse mechanisms, a number of internal routines, and the auditory nerve model. The tests should pass for each component. Failures may indicate incorrect installation or incorrect function within individual components.

Additional tests are included in the main directory. For example, test\_mechanisms.py runs a voltage clamp protocol on a selected mechanism and displays the result. test\_cells.py can run protocols on selected cell models.

Building Models
---------------

It is easiest to set up a model based on this package by using a separate directory, with links (ln -s) to CNModel and neuronvis (and possibly the mechanisms), rather than installing those as libraries in the Python environment. This helps to separate the detailed specifications from the hypotheses that are being examined in the modeling.



