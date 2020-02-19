Changes
=======

01 May 2019, 04 June 2019

This version of cnmodel runs under Python3.6 or Python3.7, using Neuron 7.6. New features include a method for changing the data tables on the fly without editing the original tables, and a tool for fitting Exp2Syn "simple" PSCs to the multisite PSC data (or, potentially, to experimental data) to get parameters for the synapse description table.

The code base has been modified for Python 3. Functionally, the main internal change is that the parameters for the cells are (almost) completely removed to the data tables. All tests currently pass, but in a few cases are very close but not identical to the original Python 2.7 version (aka branch "master-Python27"). The source of one set of discrepancies has been traced to an error in a .mod file (a variable was declared in both the RANGE and GLOBAL lists).


About CNModel
=============

CNModel is a Python-based interface to NEURON models of cochlear nucleus neurons, synapses, network connectivity. To drive the model with sound stimuli, the Zilany et al (2010, 2014) auditory periphery model is used to generate auditory nerve spike trains (either via the "cochlea" Python package or by the original MATLAB model; see below). The overall goal is to provide a platform for modeling networks of cochlear nucleus neurons with different levels of biological realism in the context of an accurate simulation of auditory nerve input.

At the lowest level are NEURON-based implementations of ion channels and synapses. Ion channel models for potassium channels are derived largely from the measurements and models of Rothman and Manis (2003abc), and Kanold and Manis (1999, 2001); other channels are derived or modified from the literature. Cell models are in turn based on the insertion of ion channels in densities based on measurements in guinea pig and mouse. The "point" somatic cell models, which form the base set of models in CNModel, replicate the data reported in the original papers. 

The postsynaptic receptor/conductance models are based on kinetic models of glutamate (Raman and Trussell, 1992) and glycinergic receptors, as adjusted to match measurements of synaptic conductances from the mouse cochlear nucleus collected in Xie and Manis, 2013. The glutamate receptor models include desensitization and the effects of internal polyamine receptor block, based on the kinetic scheme of Woodhull (1982).

The presynaptic release model includes a multisite, probabilistic synapse that includes time-dependent changes in release probability based on the Dittman, Kreitzer and Regehr (J Neurosci. 2000 Feb 15;20(4):1374-85) kinetic scheme. Although detailed, this model is computationally expensive and likely not suitable for large scale network simulations. Other simpler models of synapses are also included.

Network connectivity may be defined programmatically, or based on a table of connectivity patterns. A table with estimates derived from the literature is included in the source. 

Included in this package is a set of test suites for different components. An overriding set of unit tests is available to confirm performance of the base models against a set of reference runs, several of which are in turn directly traceable to the original manuscripts. The test suites are useful in verifying performance of the model after modifications of the code or changes in the installation (upgrades of Python or Matlab, for example). 

A manuscript describing this package has been published:
--------------------------------------------------------

    Paul B. Manis, Luke Campagnola,
    A biophysical modelling platform of the cochlear nucleus and other auditory circuits: 
    From channels to networks,
    Hearing Research,
    Volume 360,
    2018,
    Pages 76-91,
    ISSN 0378-5955,
    https://doi.org/10.1016/j.heares.2017.12.017.
    Open Access: http://www.sciencedirect.com/science/article/pii/S037859551730360X

If you use this package, we would appreciate it if you cite our work in any publications or abstracts.


Installation requirements
-------------------------
This package depends on the following:

1. Python 3.6 or 3.7 with numpy (1.14.3 or later), scipy (1.1.0 or later), lmfit (0.9.11 or later), matplotlib (3.0.3), faulthandler, and pyqtgraph (0.11.0). The cochlea module requires pandas as well. 
   An Anaconda install with the appropriate scientific packages works well::
       
       conda install python=3.6 pyqt pyqtgraph matplotlib numpy scipy pandas pytest cython
       pip install resampy
       pip install lmfit
       pip install cochlea
       
       or:
       
       conda create --name py3mpl3 python=3.6 pyqt pyqtgraph matplotlib=3 numpy scipy pandas pytest cython
       pip install resampy
       pip install lmfit
       pip install cochlea
       
      
       (Note that under MacOSX, python 3.7 is usable, including with Matlab R2019a; the Windows version of Matlab R2018b is restricted
           to python 3.6)


2. A Python-linked version of NEURON (www.neuron.yale.edu). The code has been tested with NEURON 7.5 and 7.6. We recommend
getting the most recent version of NEURON and recompiling the .mod files in the mechanisms directory.

3. A C compiler (gcc). Needed for compilation of mechanisms for NEURON.

4. The Zilany et al (JASA 2014) auditory periphery model.

This can be provided one of two ways:
    
   * The Python-based cochlea model by Rudnicki and Hemmert at https://github.com/mrkrd/cochlea. 
     This is probably best installed via pip per the instructions on the PyPi site: ``pip install cochlea``.
   * The original MATLAB-based Zilany model; requires MATLAB 2011 or later. A C compiler will also
     be needed to build this model. The model should be placed in the directory 
     ``cnmodel/cnmodel/an_model/model``, with the following components: ANmodel.m, complex.c, complex.h, 
     complex.hpp, ffGn.m, fituaudiogram2.m, mexANmodel.m, model_IHC.c, model_Synapse.c, 
     and the THRESHOLD_ALL_* files. When cnmodel attempts to access this code the first time, 
     it will perform the necessary compilation.
   
5. neuronvis (optional) available at https://github.com/campagnola/neuronvis or (a newer version) https://github.com/pbmanis/neuronvis).
   This provides 3D visualization for morphology, and is independent of cnmodel. neuronvis will require: mayavi, matplotlib, and pyqtgraph.

Once CNModel has been downloaded, go to the directory, and make sure that you are using the right branch ("Python3")::
        
        $ cd cnmodel
        $ git branch           # the active branch will have "*" next to it
        $ git checkout Python3 #(if necessary)

After the code is installed, enter the cnmodel directory and compile the NEURON mod files::

        $ nrnivmodl cnmodel/mechanisms

This will create a directory ("x86_64" or "special") in the top cnmodel directory with the compiled mechanisms.

    Under Windows 10, use::

         $ mknrndll cnmodel\mechanisms

to do the same thing. 

    Finally, go into the cnmodel directory and run::
    
        python setup.py develop
        or:
        python setup.py install

We prefer the "develop" method, as it allows you to change the code in the cnmodel directory if necessary, without re-running the setup command.


Windows Notes:
--------------

1. For more detailed information on setup in a Windows environment for Python 2.7, see the file Windows_setup.md. Thanks to Laurel Carney for prompting the generation of this set of instructions, and for identifying issues on Windows. A similar setup should work for Python 3.6+.

2. Manually compile the mex files for the Zilany et al model. In Matlab, go to the an_model/models folder, and use mexANmodel.m to compile the files. Then, add the an_model/model folder to the Matlab path, so that it can find the files when needed.

3. Under Windows, it may be best to use the standard Windows command terminal rather than the "bash" terminal provided by NEURON, at least to run the Python scripts.


Testing
-------

Make sure you are in the cnmodel directory, and that you have selected the right environment in Anaconda (in 
my case, this is usually an environment called py3mpl3 - python 3 with matplotlib 3).

After the code is installed, enter the cnmodel directory and compile the NEURON mod::

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


Usage
-----
CNModel is meant to be used as an imported package under Python. See the files in the examples directory to see how this is done. Typically, we create a separate directory (a "simulation" directory) that holds the code that uses cnmodel for simulations, at the same level as cnmodel or elsewhere (do not place the simulation directory inside cnmodel).

The data tables in the cnmodel/data directory (synapses, ionchannels, populations, connectivity) should not be modified. If it is desired to change the parameters specified in these tables, it is best to copy them into the "simulation" directory in a separate path, and modify them there. The data tables can then be used as follows::

        from cnmodel import data
        import data_XM13nacncoop as CHAN  # where data_XM13nacncoop.py is a modified table in the simulation directory
        # The following takes the table named "XM13nacncoop_channels" from the CHAN table,
        # and overwrites the original table "XM13nacncoop_channels" that is in ionchannels.py. The original file in cnmodel is
        # not modified, only the data in memory. 
        changes = data.add_table_data('XM13nacncoop_channels', row_key='field', col_key='model_type',
                       species='mouse', data=CHAN.ChannelData)
        # The following takes the table indicating how the channel compartments should be decorated from the ChannelComparments
        # table, overwriting the original named table in ionchannels.py
        changes_c = data.add_table_data('XM13nacncoop_channels_compartments', row_key='parameter', col_key='compartment',
                species='mouse', model_type='II', data=CHAN.ChannelCompartments)
        # now print out what was changed!
        data.report_changes(changes)
        data.report_changes(changes_c)

That is all that it takes. Note the following: 

1. There are some limitations as to which parameters can be changed. as some paramaters, such as rate constants for the receptors and ion channels, are specified in the .mod files and are not exposed externally. 

2. The connectivity data table can be modified to represent a particular pattern of connectivity, and the populations data table can be modified to change the relative numbers of cells.
        
3. The data tables are very strict about column alignment. The first character of the column title and the each of the values in that column must line up directly. It is best/easiest to edit these tables in a programming editor with fixed width fonts and the ability to perform column-based insertions. Changes to the data tables should be annotated appropriately.

4. Channels and receptors are all specified as NEURON .mod files. Adding new mechanisms to a cell will require modification of the code to recognize the mechanisms at several points. It is especially to handle this in cnmodel/cells.py, where knowledge of channel names is needed to compute initial states; in the cell code itself where the channels are actually inserted, and in the relevant data tables. Specific naming conventions should be followed to simplify integration. Contact the authors for help.

Adding new cell types
---------------------

To add a new cell type, it is necessary to:
    
1. Create a source file in cnmodel/cells, likely based on the bushy.py source, renaming variables as necessary. The main routines in the class however, should maintain their present names and calling parameters.
    
1. Add the values for the cells to the data tables (all tables will need to be updated with new columns for the cell type).

1. Run the model and make sure the new cell type is performing as desired. Target parameters should be identified and verified against the model.

1.  Update the unit tests to include the new cell type.





MATLAB (R)
------
This version has been tested with the MATLAB AN model of Zilany et al., 2014. 
Before using, you will need to compile the C code in an_model using Matlab's mex tool. First however, it *may* be necessary to change the following code:

In model_Synapse.c (cnmodel/an_model/model):

Change (line 63 in the source)::

	$ int    nrep, pxbins, lp,  outsize[2], totalstim;

to::

	$ int    nrep, pxbins, lp,  totalstim;
    $ size_t outsize[2];
    
Likewise, in model_IHC.c, change::

	$ int    nrep, pxbins, lp,  outsize[2], totalstim, species;

to::

	$ int    nrep, pxbins, lp,  totalstim, species;
    $ size_t outsize[2];

Then, in Matlab, go to the cnmodel/an_model/model directory, and run::

    $ mexANmodel

Then, cd to an_model and run::
    
    $ testANmodel    
    
to confirm that the model is installed and working.
(You may need to add the model directory to the Matlab path.)


Figures
-------

The data for most of the figures in the manuscript (Manis and Campagnola, Hearing Research 2018) can be generated using the bash script "figures.sh" in the examples subdirectory. 
From the main cnmodel directory::

    $ ./examples figures.sh fignum

where fignum is one of 2a, 2b, 2c, 3, 4, 5, 6a, 6b, or 7.

Note that Figure 7 may take several **hours** to generate.

Example code and tests
----------------------

A number of additional tests are included in the examples directory.

    
- `test_an_model.py` verifies that the auditory nerve model can be run. If necessary, it will compile (using MEX) the mechanisms for matlab. 
- `test_ccstim.py` tests the generation of different stimulus waveforms by the pulse generator module.
- `test_cells.py` runs different cell models in current or voltage clamp. 
  
  Usage::
      
      test_cells.py celltype species[-h] [--type TYPE] [--temp TEMP] [-m MORPHOLOGY]
                    [--nav NAV] [--ttx] [-p PULSETYPE] [--vc | --cc | --rmp]
                    
  For example: ``python test_cells.py bushy mouse --cc --temp 34``

                  
- `test_cells.py` can run protocols on selected cell models.
  Usage:: 
    
        test_cells.py [-h] [--type TYPE] [--temp TEMP] [-m MORPHOLOGY]
                      [--nav NAV] [--ttx] [-p PULSETYPE] [--vc | --cc | --rmp]
                      celltype species

- `test_circuit.py` tests the generation of circuits with populations of cells. No simulations are run.
- `test_decorator.py` generates an IV curve for the reconstructed cell LC_bushy.hoc (Figure 5B,C)
- `test_mechanisms.py` runs a voltage clamp I/V protocol on a selected mechanism and displays the result.
  Usage:: 
       
         python test_mechanisms.py <mechname>
           
  Available channel mechanisms:
              
   ========== ========= ========== ============= ==================
    CaPCalyx   KIR       bkpkj      hcno          hcnobo           
    hh         hpkj      ihpyr      ihsgcApical   ihsgcBasalMiddle 
    ihvcn      jsrna     k_ion      ka            kcnq             
    kdpyr      kht       kif        kis           klt              
    kpkj       kpkj2     kpkjslow   kpksk         leak             
    lkpkj      na        naRsg      na_ion        nacn             
    nacncoop   nap       napyr      nav11                          
   ========== ========= ========== ============= ==================

- `test_mso_inputs.py` runs a circuit that creates a point MSO neuron, innervated by bushy cells from independent "ears". This demonstrates how to construct a binaural circuit using CNModel.
- `test_physiology.py` runs a large VCN circuit that converges onto a single bushy cell. This run can take a long time. The output was used to create Figure 7 of the manuscript.
- `test_populations.py` tests synaptic connections between two cell types. Usage::
    
      python test_populations.py <pre_celltype> <post_celltype>
      
- `test_sgc_input_phaselocking.py` tests phase locking with SGC inputs to a bushy cell.
- `test_sgc_input_PSTH.py` shows SGC inputs and postsynaptic bushy cell PSTHs.
- `test_sgc_input.py` demonstrates SGC input to a VCN bushy cell.
- `test_simple_synapses.py` tests simple Exp2Syn inputs to different cell types. Usage::
    
      python test_synapses.py <pre_celltype> <post_celltype>
      
  Supported cell types: sgc, bushy, tstellate, dstellate, tuberculoventral, pyramidal
- `test_sound_stim.py` generates spike trains from the selected model (cochlea, matlab) and plots rate-intensity functions for the 3 different SR groups.
- `test_sounds.py` generates waveforms for different kinds of sounds included in the sounds class.
- `test_synapses.py` evokes spikes in a presynaptic cell while recording the postsynaptic potential. Usage::
    
      python test_synapses.py <pre_celltype> <post_celltype>
      
  Supported cell types: sgc, bushy, tstellate, dstellate
- `toy_model.py` generates IV plots for each of the principal point cell types included in CNModel. This is the code that generates Figure 3 of the manuscript.

Potential Issues and Solutions
------------------------------

1.  Occasionally one of the AN spike train files, which are stored in the directory `cnmodel/an_model/cache`, become locked. This can occur if the calling routines (e.g., simulation runs) are aborted (^C, ^Z) in the middle of a transaction accessing the cache file, or perhaps during when parallel processing is enabled and a routine fails or is aborted. In this case, a file with the extension ``".lock"`` exists, which prevents the an_model code from accessing the file. The ``".lock"`` file needs to be deleted from the cache directory. Because the cache directory contains an hierarchical arrangement of subdirectories, and can be populated with thousands of files after a few runs requiring many auditory nerve datasets, finding the lock file can be somewhat tedious. The following should help under Unix:
    
  *  First, print a list of the locked files::
      
          - find /path/to/cache -name '*.lock'
    
  * Where /path/to/cache may be something like `cnmodel/an_model/cache`. 
    There is most likely only one such file in the diretory.

  * Next, to delete the files::
  
      - find /path/to/cache -name '*.lock' -delete
       
  * Under Windows (and other OS's), you should be able do accomplish the same thing
    with the File Explorer/Finder, limiting the files by extension.
    
  * An alternative (for any OS) is to take advantage of Python's pathlib module. The glob search is 
    remarkably fast (on my system, it takes under a minute to search through more than 3.5 million
    cached AN spike trains)::
    
            >>python
            > from pathlib import Path
            > gl = Path('.').rglob('*.lock')
            > locks = list(gl) # (could do this in the next line)
            > # print(locks)  # see the lock files
            > for g in locks:  # now remove the lock files
            >    g.unlink()
            >
   
References
----------

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

