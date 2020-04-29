About CNModel
=============

CNModel is a Python-based interface to NEURON models of cochlear nucleus neurons, synapses, network connectivity. To drive the model with sound stimuli, the Zilany et al (2010, 2014) auditory periphery model is used to generate auditory nerve spike trains (either via the "cochlea" Python package or by the original MATLAB model; see below). The overall goal is to provide a platform for modeling networks of cochlear nucleus neurons with different levels of biological realism in the context of an accurate simulation of auditory nerve input.

At the lowest level are NEURON-based implementations of ion channels and synapses. Ion channel models for potassium channels are derived largely from the measurements and models of Rothman and Manis (2003abc), and Kanold and Manis (1999, 2001); other channels are derived or modified from the literature. Cell models are in turn based on the insertion of ion channels in densities based on measurements in guinea pig and mouse. The "point" somatic cell models, which form the base set of models in CNModel, replicate the data reported in the original papers. 

The postsynaptic receptor/conductance models are based on kinetic models of glutamate (Raman and Trussell, 1992) and glycinergic receptors, as adjusted to match measurements of synaptic conductances from the mouse cochlear nucleus collected in Xie and Manis, 2013. The glutamate receptor models include desensitization and the effects of internal polyamine receptor block, based on the kinetic scheme of Woodhull (1982).

The presynaptic release model includes a multisite, probabilistic synapse that includes time-dependent changes in release probability based on the Dittman, Kreitzer and Regehr (J Neurosci. 2000 Feb 15;20(4):1374-85) kinetic scheme. Although detailed, this model is computationally expensive and likely not suitable for large scale network simulations. Other simpler models of synapses are also included.

Network connectivity may be defined programmatically, or based on a table of connectivity patterns. A table with estimates derived from the literature is included in the source. 

Included in this package is a set of test suites for different components. An overriding set of unit tests is available to confirm performance of the base models against a set of reference runs, several of which are in turn directly traceable to the original manuscripts. The test suites are useful in verifying performance of the model after modifications of the code or changes in the installation (upgrades of Python or Matlab, for example).


Changes
-------

February 19, 2020:

The Python 3 version of CNModel has been made the default version.


01 May 2019, 04 June 2019:

This version of cnmodel runs under Python3.6 or Python3.7, using Neuron 7.6. New features include a method for changing the data tables on the fly without editing the original tables, and a tool for fitting Exp2Syn "simple" PSCs to the multisite PSC data (or, potentially, to experimental data) to get parameters for the synapse description table.

The code base has been modified for Python 3. Functionally, the main internal change is that the parameters for the cells are (almost) completely removed to the data tables. All tests currently pass, but in a few cases are very close but not identical to the original Python 2.7 version (aka branch "master-Python27"). The source of one set of discrepancies has been traced to an error in a .mod file (a variable was declared in both the RANGE and GLOBAL lists).

