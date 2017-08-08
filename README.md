About CNModel
=============

CNModel is a Python-based interface to NEURON models of cochlear nucleus neurons, synapses, network connectivity. To drive the model with sound stimuli, the Zilaney et al (2010, 2014) auditory periphery model is used to generate auditory nerve spike trains (either via the "cochlea" Python package or by the original MATLAB model; see below). The overall goal is to provide a platform for modeling networks of cochlear nucleus neurons with different levels of biological realism in the context of an accurate simulation of auditory nerve input.

At the lowest level are NEURON-based implementations of ion channels and synapses. Ion channel models for potassium channels are derived largely from the measurements and models of Rothman and Manis (2003abc), and Kanold and Manis (1999, 2001); other channels are derived or modified from the literature. Cell models are in turn based on the insertion of ion channels in densities based on measurements in guinea pig and mouse. The "point" somatic cell models, which form the base set of models in CNModel, replicate the data reported in the original papers. 

The postsynaptic receptor/conductance models are based on kinetic models of glutamate (Raman and Trussell, 1992) and glycinergic receptors, as adjusted to match measurements of synaptic conductances from the mouse cochlear nucleus collected in Xie and Manis, 2013. The glutamate receptor models include desensitization and the effects of internal polyamine receptor block, based on the kinetic scheme of Woodhull (1982).

The presynaptic release model includes a multisite, probabilistic synapse that includes time-dependent changes in release probability based on the Dittman, Kreitzer and Regehr (J Neurosci. 2000 Feb 15;20(4):1374-85) kinetic scheme. Although detailed, this model is computationally expensive and likely not suitable for large scale network simulations. Other simpler models of synapses are also included.

Network connectivity may be defined programmatically, or based on a table of connectivity patterns. A table with estimates derived from the literature is included in the source. 

Included in this package is a set of test suites for different components. An overriding set of unit tests is available to confirm performance of the base models against a set of reference runs, several of which are in turn directly traceable to the original manuscripts. The test suites are useful in verifying performance of the model after modifications of the code or changes in the installation (upgrades of Python or Matlab, for example). 

Installation requirements
-------------------------
This package depends on the following:

   1. Python 2.7.10 with numpy (1.11), scipy (0.19.0), lmfit (0.9.6), pyqt4 (4.11.4), and pyqtgraph (0.9.10). 
      An Anaconda install with the appropriate scientific packages works well. lmfit is best obtained via pip
      to install the latest versions. This package is not yet compatible with Python 3.x.
   2. A Python-linked version of NEURON (www.neuron.yale.edu). The code has been tested with NEURON 7.3 and 7.4.
   3. A C compiler (gcc). Needed for compilation of mechanisms for NEURON.
   4. The Zilany et al (JASA 2014) auditory periphery model. This can be provided one of two ways:
      1. The original MATLAB-based Zilany model; requires MATLAB 2011 or later. A C compiler will also be needed to build this model.
      2. The Python-based cochlea model by Rudnicki and Hemmert (https://github.com/mrkrd/cochlea; can be installed via pip)
   5. neuronvis (optional; available at https://github.com/campagnola/neuronvis or https://github.com/pbmanis/neuronvis).
      This provides 3D visualization for morphology.


Testing
-------

Before testing, enter the cnmodel directory and compile the NEURON mod files:

    $ nrnivmodl cnmodel/mechanisms

This will create a directory ("x86_64" or "special") in the top cnmodel directory. At that point

    $ python examples/toy_model.py
     
should generate a plot with several sets of traces showing responses of individual neuron models to depolarizing and hyperpolarizing current steps.

The test suite can be run as:

    $ python test.py

This will test each of the models against reference data, the synapse mechanisms, a number of internal routines, and the auditory nerve model. The tests should pass for each component. Failures may indicate incorrect installation or incorrect function within individual components.

Additional tests are included in the examples directory. For example:
    
    * `test_mechanisms.py` runs a voltage clamp I/V protocol on a selected mechanism and displays the result.
    * `test_cells.py` can run protocols on selected cell models.
    * `test_synapses.py` evokes spikes in a presynaptic cell while recording the postsynaptic potential.
    
    
References:
-----------

1.   Cao XJ, Oertel D. The magnitudes of hyperpolarization-activated and
low-voltage-activated potassium currents co-vary in neurons of the ventral
cochlear nucleus. J Neurophysiol. 2011 Aug;106(2):630-40. doi:
10.1152/jn.00015.2010. Epub 2011 May 11. PubMed PMID: 21562186; PubMed Central
PMCID: PMC3154804.

2.   Cao XJ, Oertel D. Auditory nerve fibers excite targets through synapses that
vary in convergence, strength, and short-term plasticity. J Neurophysiol. 2010
Nov;104(5):2308-20. doi: 10.1152/jn.00451.2010. Epub 2010 Aug 25. PubMed PMID:
20739600; PubMed Central PMCID: PMC3350034.

3.   Dittman JS, Kreitzer AC, Regehr WG. Interplay between facilitation, depression,
and residual calcium at three presynaptic terminals. J Neurosci. 2000 
Feb 15;20(4):1374-85. PubMed PMID: 10662828.

1. Isaacson JS, Walmsley B. Counting quanta: direct measurements of transmitter
release at a central synapse. Neuron. 1995 Oct;15(4):875-84.

4. Kanold PO, Manis PB. A physiologically based model of discharge pattern
regulation by transient K+ currents in cochlear nucleus pyramidal cells. J
Neurophysiol. 2001 Feb;85(2):523-38. PubMed PMID: 11160490.

5.   Kanold PO, Manis PB. Transient potassium currents regulate the discharge
patterns of dorsal cochlear nucleus pyramidal cells. J Neurosci. 1999 Mar
15;19(6):2195-208. PubMed PMID: 10066273.

6.   Liu Q, Manis PB, Davis RL. Ih and HCN channels in murine spiral ganglion
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

