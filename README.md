cnmodel
=======

cnmodel is a Python-based interface to NEURON models of cochlear nucleus neurons, synapses, network connectivity, and the Zilany et al (2009, 2010, 2014) model of auditory nerve fiber spike trains to acoustic stimuli. The overall goal is to provide a platform for modeling networks of cochlear nucleus neurons with different levels of biological realism in the context of an accurate simulation of auditory nerve input.

At the lowest level are NEURON-based implementations of ion channels and synapses. Ion channel models for potassium channels are derived largely from the measurements and models of Rothman and Manis (1993, 2003abc), and Kanold and Manis (1999, 2001); other channels are derived or modified from the literature (for example, Liu et al, 2014). Cell models are in turn based on the insertion of ion channels in densities based on measurements in guinea pig (Rothman and Manis, 2003b) and mouse (Cao and Oertel, 2010, 2011). The "point" somatic cell models for VCN neurons (bushy, t-stellate), which form the base set of models in cnmodel, replicate the data reported in the original papers of Rothman and Manis (b,c). 

The postsynaptic receptor/conductance models are based on kinetic models of glutamate (Raman and Trussell, 1992) and glycinergic receptors, as adjusted to match measurements of synaptic conductances from the mouse cochlear nucleus collected in Xie and Manis, (2013). The glutamate receptor models include desensitization and the effects of internal polyamine receptor block, based on the kinetic scheme of Woodhull (1973).

The presynaptic release model includes a multisite, probabilistic synapse that includes time-dependent changes in release probability based on the Dittman, Kreiger and Regehr (2000) kinetic scheme, as well as variable quantal latency (based on Isaacson and Walmsely, 1995) and time-dependent latency (based on unpublished measurements). Because it is detailed, this model is computationally expensive and not suitable for large scale network simulations. Aspects of the release model can be "turned off" if they are not deeemed important, to speed up simulations.

Network connectivity may be defined programatically, or based on a table of connetivity patterns. A table with estimates and guesses derived from the literature is included in the source. 

Included in this package is a set of test suites for different components. An overriding set of unit tests is available to confirm performance of the base models against a set of reference runs, several of which are in turn directly traceable to the original manuscripts. The test suites are useful in verifying performance of the model after modifications of the code or changes in the installation (upgrades of Python or Matlab, for example). 

Installation requirements
-------------------------
This package depends on the following:

   1. Python 2.7.10 with at least numpy, scipy, Qt4.8 and pyqtgraph. An Anaconda install with the appropriate scientific packages works well. This project has not been tested in Python 3.x.
   2. A Python-linked version of NEURON (www.neuron.yale.edu). The code has been tested with NEURON 7.3 and 7.4.
   3. A C compiler (gcc). Needed for compilation of C files for matlab, and compilation of the mechanisms for NEURON.
   4. Matlab (2011 or later). Matlab is needed for the auditory nerve model to run. The C compiler is needed to compile the C files in the auditory nerve model. 
   5. neuronvis (available at https://github.com/campagnola/neuronvis). This provides a library that can read and visualize hoc files.
   6. lmfit (avaliable at https://lmfit.github.io/lmfit-py/index.html). You need to use a version < 0.9.0 (the code works with 0.8.0); API changes in lmfit for 0.9.0 break the handling of the fit parameters. This is not fixed yet. Install with "pip install 'lmfit==0.8.0'" and check that the correct version is imported when running Python (import lmfit; print lmfit.__version__). 

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

Current Issues with test suites
-------------------------------

Currently, the tests indicate multiple failiures (many in master branch). Most of these have been fixed in the "hoc-incorporation" branch (4 failures). These 4 failures represent two sets of failures, and each set indicates one real failure and one cascading error in the test suite itself. The first failure (cnmodel/cells/tests/test_cells.py::test_octopus FAILED; cnmodel/cells/tests/test_cells.py::test_pyramidal FAILED) means that the archived test data for the octopus cell does not precisely meet the results of current runs. The model is qualitatively the same, but quantitatively differnt. This is being investigated. The pyramidal cell test failure is an bug in the test suite that follows from the previously failed test.

The same issue exists for the two synapse tests that fail: (cnmodel/synapses/tests/test_synapses.py::test_dstellate_bushy FAILED
cnmodel/synapses/tests/test_synapses.py::test_dstellate_tstellate FAILED). The first is a quantative difference that has not yet been investigated (the results *look* right, but are quantitatively not matching to the archived reference run); the second is a cascadeing error in the test suite from the first. 


References:

1.   Cao XJ, Oertel D. The magnitudes of hyperpolarization-activated and
low-voltage-activated potassium currents co-vary in neurons of the ventral
cochlear nucleus. J Neurophysiol. 2011 Aug;106(2):630-40. doi:
10.1152/jn.00015.2010. Epub 2011 May 11. PubMed PMID: 21562186; PubMed Central
PMCID: PMC3154804.

2.   Cao XJ, Oertel D. Auditory nerve fibers excite targets through synapses that
vary in convergence, strength, and short-term plasticity. J Neurophysiol. 2010
Nov;104(5):2308-20. doi: 10.1152/jn.00451.2010. Epub 2010 Aug 25. PubMed PMID:
20739600; PubMed Central PMCID: PMC3350034.

3.   Dittman JS, Kreitzer AC, Regehr WG. Interplay between facilitation, depression, and residual calcium at three presynaptic terminals. J Neurosci. 2000 Feb 15;20(4):1374-85. PubMed PMID: 10662828.

1. Isaacson JS, Walmsley B. Counting quanta: direct measurements of transmitter
release at a central synapse. Neuron. 1995 Oct;15(4):875-84. PubMed PMID:
7576636.

4.   Kanold PO, Manis PB. A physiologically based model of discharge pattern
regulation by transient K+ currents in cochlear nucleus pyramidal cells. J
Neurophysiol. 2001 Feb;85(2):523-38. PubMed PMID: 11160490.


5.  Kanold PO, Manis PB. Transient potassium currents regulate the discharge
patterns of dorsal cochlear nucleus pyramidal cells. J Neurosci. 1999 Mar
15;19(6):2195-208. PubMed PMID: 10066273.

6.   Liu Q, Manis PB, Davis RL. I h and HCN channels in murine spiral ganglion
neurons: tonotopic variation, local heterogeneity, and kinetic model. J Assoc Res
Otolaryngol. 2014 Aug;15(4):585-99. doi: 10.1007/s10162-014-0446-z. Epub 2014 Feb
21. Erratum in: J Assoc Res Otolaryngol. 2014 Aug;15(4):601. PubMed PMID:
24558054; PubMed Central PMCID: PMC4141436.

7.   Raman IM, Trussell LO. The kinetics of the response to glutamate and kainate
in neurons of the avian cochlear nucleus. Neuron. 1992 Jul;9(1):173-86. PubMed
PMID: 1352983.

8.   Rothman JS, Manis PB. The roles potassium currents play in regulating the
electrical activity of ventral cochlear nucleus neurons. J Neurophysiol. 2003
Jun;89(6):3097-113. PubMed PMID: 12783953.

9.  Rothman JS, Manis PB. Kinetic analyses of three distinct potassium
conductances in ventral cochlear nucleus neurons. J Neurophysiol. 2003
Jun;89(6):3083-96. PubMed PMID: 12783952.

10.   Rothman JS, Manis PB. Differential expression of three distinct potassium
currents in the ventral cochlear nucleus. J Neurophysiol. 2003 Jun;89(6):3070-82.
PubMed PMID: 12783951.


11.   Rothman JS, Young ED, Manis PB. Convergence of auditory nerve fibers onto
bushy cells in the ventral cochlear nucleus: implications of a computational
model. J Neurophysiol. 1993 Dec;70(6):2562-83. PubMed PMID: 8120599.

12.   Woodhull AM. Ionic blockage of sodium channels in nerve. J Gen Physiol. 1973
Jun;61(6):687-708. PubMed PMID: 4541078; PubMed Central PMCID: PMC2203489.

13.   Xie R, Manis PB. Target-specific IPSC kinetics promote temporal processing in 
auditory parallel pathways. J Neurosci. 2013 Jan 23;33(4):1598-614. doi:
10.1523/JNEUROSCI.2541-12.2013. PubMed PMID: 23345233; PubMed Central PMCID:
PMC3737999.

14.   Zilany MS, Bruce IC, Carney LH. Updated parameters and expanded simulation
options for a model of the auditory periphery. J Acoust Soc Am. 2014
Jan;135(1):283-6. doi: 10.1121/1.4837815. PubMed PMID: 24437768; PubMed Central
PMCID: PMC3985897.

15.   Zilany MS, Carney LH. Power-law dynamics in an auditory-nerve model can
account for neural adaptation to sound-level statistics. J Neurosci. 2010 Aug
4;30(31):10380-90. doi: 10.1523/JNEUROSCI.0647-10.2010. PubMed PMID: 20685981;
PubMed Central PMCID: PMC2935089.


16.   Zilany MS, Bruce IC, Nelson PC, Carney LH. A phenomenological model of the
synapse between the inner hair cell and auditory nerve: long-term adaptation with
power-law dynamics. J Acoust Soc Am. 2009 Nov;126(5):2390-412. doi:
10.1121/1.3238250. PubMed PMID: 19894822; PubMed Central PMCID: PMC2787068.





