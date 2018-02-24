Windows setup
=============

(as of 24 Feb 2018, tested on Windows 10, )
On a bare system, you will need to install the following packages:

    1. Anaconda Python 2.7 (this may also include the Microsoft Visual Studio for Python)
    2. git
    3. msvc2.7forpython (required for cochlea)
    4. external python package from PyPi: cochlea
    5. NEURON7.5 (this must be installed last)
    6. cnmodel, cloned from the main repository with git
    7. build the .dll file for the Neuron mechanisms used in cnmodel

Some installs are accomplished via a graphical interface (anaconda, git, msvc, neuron). 
Other parts of the install are done from the terminal (Winows "Command Prompt" - look in the Windows Accessories)

Follow the instructions below to do this installation. The order only partially matters - msvcforpython must be installed before cochlea can be built; python should be installed before NEURON; and git is required get the cnmodel repository.

Step 1: Install Anaconda python for Python 2.7 for your system. 
Note that NEURON is not yet compatible with Python 3.x

https://repo.continuum.io/archive/
The most recent versions of the Anaconda installer (5.x) are OK to use.
Check the box to set DOS path variables so can run from command prompt.
Install only for yourself.

If you want, create an environment for python 2.7, Using the standard windows terminal ("Command Prompt"), install the additional required pacakges::

    conda create --name models python=2.7 pyqt pyqtgraph matplotlib numpy scipy pytest faulthandler

Note: this may take a while, and will install many other packages that are required dependencies.

Now, activate the environment::

    conda activate models

Or else, just install everything to the root environment:
conda install pyqt pyqtgraph matplotlib numpy scipy pytest faulthandler

If you are working with a previous Anaconda installation, you might need to do these steps::

    conda update --all
    conda update qt

The following must be retrieved from the PyPi repository (I had problems with the one in the anaconda repo)::

    pip install lmfit

Step 2: Install git (a widely used source control module)::

    https://git-scm.com/downloads

Step 3: Install MS visual studio for python::

    https://www.visualstudio.com/vs/python/

Step 4: Install cochlea (Rudnicki and Hemmert's implementation of several models under python)::

    pip install cochlea

This should build quickly and have no errors. 

Step 5: Install NEURON7.5 for mswindows from::
    www.neuron.yale.edu. 

Restart the command window so that the paths to the python installation are updated.
Make sure that python and neuron are installed correctly together by testing the ability to import NEURON::

    (base) pbmanis: Python $ python
    Python 2.7.14 |Anaconda custom (64-bit)| (default, Dec  7 2017, 11:07:58)
    [GCC 4.2.1 Compatible Clang 4.0.1 (tags/RELEASE_401/final)] on darwin
    Type "help", "copyright", "credits" or "license" for more information.
    >>> from neuron import h
    NEURON -- VERSION 7.5 master (6b4c19f) 2017-09-25
    Duke, Yale, and the BlueBrain Project -- Copyright 1984-2016
    See http://neuron.yale.edu/neuron/credits

    >>>
    
(^Z to exit)

Step 6: Go to where you want to put the model code repository, and clone the repo::

    git clone https://github.com/cnmodel/cnmodel.git

Jump into the cnmodel directory.
Currently, you will need to switch the active branch to Qt5::

    git checkout Qt5

(I hope to remove this step soon).

Step 7: Find the mknrndll gui from where you installed neuron with the Explorer (the program in the nrn folder). Run the program, and select the directory cnmodel\cnmodel\mechanisms, then build.
This will make nrnmech.dll in the mechanisms directory.
Copy the nrnmech.dll into the main cnmodel directory.

For Windows, you may also need to copy the example files out of the examples directory into the main cnmodel directory (this is where they lived before we cleaned up the organization, and probably never tested against windows again).

You can run "python setup.py build install" in the main cnmodel directory to make a library version. This will make cnmodel accessible as an import into python from any location. The only catch is that the nrnmech.dll must be placed locally in the directory that you are running from, and it may lead to problems in the future.


Running
=======

The following should be all that is needed once everything is set up.

Go the the cnmodel main directory.

if you have an environment, activate it to make sure in the right anaconda python environment

Now, these scripts should run just fine::

    python toy_model.py 
    python test_mechanisms.py klt
    python test_synapses.py sgc bushy
    python test_phaselocking.py   will fail as it needs the matlab AN model.

You should also be able to run the full battery of tests::

    python test.py

All the tests should pass (except maybe dstellate->dstellate; this passes under MacOSX so I'll have to investigate why it fails under windows). MacOSX test result is listed below.

Notes
=====
1. You cannot use the "bash" terminal window that comes with neuron - it doesn't set the paths correctly for access to the anaconda python.

2. Some of the graphical displays are not correctly sized in Windows. You may need to expand the window to see all the plots (specifically, toy_model has this problem).

3. You may encounter failures to load PyQt4 (such as, from PyQt4 import QtGui). If so, comment out the import line. Find where QtGui (or QtCore) is used and replace it with pg.QtGui (or pg.QtCore). This lets pyqtgraph resolve the appropriate Qt library. Report any such occurances to me (pmanis@med.unc.edu)

4. The test suite (python test.py) may fail on one test: dstellate -> dstellate synapses (as of 23Feb2018), and if you do not have matlab, the matlab test will be skipped.

5. Developing programs to use cnmodel: Under mac os x and linux, I use a separate directory for the code, and create a symbolic link to the top cnmodel directory for the imports. This lets me tinker with cnmodel while developing external code. However, I have not tried this under windows. Using the setup.py/install is the next best thing, although any changes to cnmodel must be followed by a setup/install call to update the package.



Things solved with Windows
==========================

1. The test suite works and most tests pass (python test.py).
2. Qt5 and current anaconda install are acceptable (no longer trying to force Qt4).


Potential problems
==================
If a run fails because it cannot find pylibrary for an import, you can get pylibrary as follows. However, the only place that pylibrary appears to be used is in a routine that is not currently used by cnmodel.

Go to a reasonable directory (not in cnmodel) - I use Desktop\Python for all of these then get pylibrary::

    git clone https://github.com/pbmanis/pylibrary
    cd pylibrary
    python setup.py build install  # puts pylibrary in the anaconda environment models only.
    cd ..
    
Pylibrary is only needed for a few of the example files.



Test results
============

Mac OSX (24 Feb 2018)

::
    (base) pbmanis: cnmodel [qt5+]$ python test.py
    NEURON -- VERSION 7.5 master (6b4c19f) 2017-09-25
    Duke, Yale, and the BlueBrain Project -- Copyright 1984-2016
    See http://neuron.yale.edu/neuron/credits

    loading membrane mechanisms from x86_64/.libs/libnrnmech.so
    Additional mechanisms from files
     cnmodel/mechanisms//CaPCalyx.mod cnmodel/mechanisms//Gly5GC.mod cnmodel/mechanisms//Gly5PL.mod 
     cnmodel/mechanisms//Gly5State.mod cnmodel/mechanisms//Gly6S.mod cnmodel/mechanisms//Iclamp2.mod 
     cnmodel/mechanisms//NMDA.mod cnmodel/mechanisms//NMDA_Kampa.mod 
     cnmodel/mechanisms//ampa_trussell.mod cnmodel/mechanisms//bkpkj.mod 
     cnmodel/mechanisms//cabpump.mod cnmodel/mechanisms//cadiff.mod cnmodel/mechanisms//cadyn.mod 
     cnmodel/mechanisms//cap.mod cnmodel/mechanisms//capmp.mod 
     cnmodel/mechanisms//capump.mod cnmodel/mechanisms//cleftXmtr.mod 
     cnmodel/mechanisms//gly.mod cnmodel/mechanisms//gly2.mod cnmodel/mechanisms//hcno.mod 
     cnmodel/mechanisms//hcno_bo.mod cnmodel/mechanisms//iStim.mod cnmodel/mechanisms//ihpkj.mod 
     cnmodel/mechanisms//ihpyr.mod cnmodel/mechanisms//ihsgc_apical.mod 
     cnmodel/mechanisms//ihsgc_basalmiddle.mod cnmodel/mechanisms//ihvcn.mod 
     cnmodel/mechanisms//inav11.mod cnmodel/mechanisms//jsrnaf.mod 
     cnmodel/mechanisms//ka.mod cnmodel/mechanisms//kcnq.mod cnmodel/mechanisms//kdpyr.mod 
     cnmodel/mechanisms//kht.mod cnmodel/mechanisms//kif.mod cnmodel/mechanisms//kir.mod 
     cnmodel/mechanisms//kis.mod cnmodel/mechanisms//klt.mod cnmodel/mechanisms//kpkj.mod 
     cnmodel/mechanisms//kpkj2.mod cnmodel/mechanisms//kpkjslow.mod cnmodel/mechanisms//kpksk.mod 
     cnmodel/mechanisms//leak.mod cnmodel/mechanisms//multisite.mod cnmodel/mechanisms//na.mod 
     cnmodel/mechanisms//nacn.mod cnmodel/mechanisms//nacncoop.mod cnmodel/mechanisms//nap.mod c
     nmodel/mechanisms//napyr.mod cnmodel/mechanisms//pkjlk.mod cnmodel/mechanisms//rsg.mod 
     cnmodel/mechanisms//vecevent.mod
    Testing with flags: -v --tb=short cnmodel/
    =============================================== test session starts ===============================================
    platform darwin -- Python 2.7.14, pytest-3.0.7, py-1.4.33, pluggy-0.4.0 -- /Users/pbmanis/anaconda/bin/python
    cachedir: .cache
    rootdir: /Users/pbmanis/Desktop/Python/cnmodel, inifile:
    collected 36 items

    cnmodel/an_model/tests/test_cache.py::test_cache PASSED
    cnmodel/an_model/tests/test_cache.py::test_parallel PASSED
    cnmodel/cells/tests/test_cells.py::test_bushy PASSED
    cnmodel/cells/tests/test_cells.py::test_bushy21 PASSED
    cnmodel/cells/tests/test_cells.py::test_bushy_mouse PASSED
    cnmodel/cells/tests/test_cells.py::test_tstellate PASSED
    cnmodel/cells/tests/test_cells.py::test_tstellate_mouse PASSED
    cnmodel/cells/tests/test_cells.py::test_tstellatet PASSED
    cnmodel/cells/tests/test_cells.py::test_dstellate PASSED
    cnmodel/cells/tests/test_cells.py::test_dstellate_mouse PASSED
    cnmodel/cells/tests/test_cells.py::test_octopus PASSED
    cnmodel/cells/tests/test_cells.py::test_pyramidal PASSED
    cnmodel/cells/tests/test_cells.py::test_tuberculoventral PASSED
    cnmodel/cells/tests/test_cells.py::test_cartwheel PASSED
    cnmodel/cells/tests/test_cells.py::test_sgc_basal_middle PASSED
    cnmodel/cells/tests/test_cells.py::test_sgc_apical PASSED
    cnmodel/data/tests/test_db.py::test_db PASSED
    cnmodel/mechanisms/tests/test_mechanisms.py::test_max_open_probability PASSED
    cnmodel/synapses/tests/test_psd.py::test_sgc_bushy_psd PASSED
    cnmodel/synapses/tests/test_psd.py::test_sgc_tstellate_psd PASSED
    cnmodel/synapses/tests/test_psd.py::test_sgc_dstellate_psd PASSED
    cnmodel/synapses/tests/test_psd.py::test_sgc_octopus_psd PASSED
    cnmodel/synapses/tests/test_synapses.py::test_sgc_bushy PASSED
    cnmodel/synapses/tests/test_synapses.py::test_sgc_tstellate PASSED
    cnmodel/synapses/tests/test_synapses.py::test_sgc_tstellate2 PASSED
    cnmodel/synapses/tests/test_synapses.py::test_sgc_dstellate PASSED
    cnmodel/synapses/tests/test_synapses.py::test_dstellate_bushy PASSED
    cnmodel/synapses/tests/test_synapses.py::test_dstellate_tstellate PASSED
    cnmodel/synapses/tests/test_synapses.py::test_dstellate_dstellate PASSED
    cnmodel/util/tests/test_expfitting.py::test_fit1 PASSED
    cnmodel/util/tests/test_expfitting.py::test_fit2 PASSED
    cnmodel/util/tests/test_matlab.py::test_matlab PASSED
    cnmodel/util/tests/test_sound.py::test_conversions PASSED
    cnmodel/util/tests/test_sound.py::test_tonepip PASSED
    cnmodel/util/tests/test_sound.py::test_noisepip PASSED
    cnmodel/util/tests/test_stim.py::test_make_pulse PASSED

    =========================================== 36 passed in 94.05 seconds ============================================
    (base) pbmanis: cnmodel [qt5+]$




