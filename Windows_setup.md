CNMODEL Setup for Windows Systems
=================================

(as of 2 March 2018, tested on Windows 10)
On a bare system, you will need to install the following packages:

    1. Anaconda Python 2.7 (this may also include the Microsoft Visual Studio for Python)
    2. git
    3. msvc2.7forpython (required for cochlea)
    4. external python package from PyPi: cochlea
    5. NEURON7.5 (this must be installed last)
    6. cnmodel, cloned from the main repository with git
    7. build the .dll file for the Neuron mechanisms used in cnmodel

Some installs are accomplished via a graphical interface (anaconda, git, msvc, neuron). 
Other parts of the install are done from the terminal (Windows "Command Prompt" - look in the Windows Accessories)

Follow the instructions below to do this installation. The order only partially matters - msvcforpython must be installed before cochlea can be built; python should be installed before NEURON; and git is required get the cnmodel repository.

Step 1: Install Anaconda python for Python 2.7 for your system. 
Note that NEURON is not yet compatible with Python 3.x

https://repo.continuum.io/archive/
The most recent versions of the Anaconda installer (5.x) are OK to use.
Check the box to set DOS path variables so can run from command prompt.
Install only for yourself.

If you want, create an environment for python 2.7, Using the standard windows terminal ("Command Prompt"), install the additional required packages::

    conda create --name models python=2.7 pyqt pyqtgraph matplotlib numpy scipy pandas pytest faulthandler

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

Step 7: Find the mknrndll gui from where you installed neuron with the Explorer (the program in the nrn folder). Run the program, and select the directory cnmodel\cnmodel\mechanisms, then build.
This will make nrnmech.dll in the mechanisms directory.
Copy the nrnmech.dll into the main cnmodel directory.

Step 8: Run::

    python setup.py develop
    
in the main cnmodel directory to make a library version in the anaconda site-packages directory. This will make cnmodel accessible as an import into python from any location. However, you will (on Windows) need to copy the nrnmech.dll file to the directory that you are starting in. (If you import cnmodel and there is not a list of mechanisms after the Neuron banner, then nrnmech.dll was not found).

Running
=======

The following should be all that is needed once everything is set up.

Go the the cnmodel main directory.

if you set up a Python environment, activate it to make sure you have the right dependencies available.

Now, these scripts should run just fine::

    python examples/test_mechanisms.py klt
    python examples/toy_model.py 
    python examples/test_synapses.py sgc bushy
    python examples/test_phaselocking.py

You should also be able to run the full battery of tests::

    python test.py

All the tests should pass (except perhaps with the exception dstellate->dstellate; this passes under MacOSX so we will need to investigate why it fails under Windows). MacOSX and Windows test results are listed below.

Notes
=====

1. You cannot use the "bash" terminal window that comes with neuron - it doesn't set the paths correctly for access to the anaconda python.

2. Some of the graphical displays are not correctly sized in Windows. You may need to expand the window to see all the plots (specifically, toy_model has this problem).

3. Report occurrences of failures to "import PyQt4" to: https://github.com/cnmodel/cnmodel/issues.

4. The test suite (python test.py) may fail on one test: dstellate -> dstellate synapses (as of 23Feb2018), and if you do not have Matlab, the Matlab test will be skipped.


Things solved with Windows
==========================

1. The test suite works and most tests pass (python test.py).
2. Qt5 and current anaconda install are acceptable (we are no longer trying to force Qt4).
3. The prior dependency on pylibrary has been removed.


Potential problems
==================

This software has been tested primarily on Mac OSX and Linux systems. Some tests have been performed on Windows. If issues are found, please report them via github as an issue (as noted above).

Test results
============

Windows 10 (24 Feb 2018; Qt5 Branch)
------------------------------------

::
    C:\Users\pbmanis\Desktop\Python\cnmodel>python test.py > testresults.txt
    NEURON -- VERSION 7.5 master (6b4c19f) 2017-09-25
    Duke, Yale, and the BlueBrain Project -- Copyright 1984-2016
    See http://neuron.yale.edu/neuron/credits

    loading membrane mechanisms from C:\Users\pbmanis\Desktop\Python\cnmodel\nrnmech.dll
    Additional mechanisms from files
     CaPCalyx.mod Gly5GC.mod Gly5PL.mod Gly5State.mod Gly6S.mod Iclamp2.mod NMDA.mod NMDA_Kampa.mod
      ampa_trussell.mod bkpkj.mod cabpump.mod cadiff.mod cadyn.mod cap.mod capmp.mod capump.mod cleftXmtr.mod
       gly.mod gly2.mod hcno.mod hcno_bo.mod iStim.mod ihpkj.mod ihpyr.mod ihsgc_apical.mod ihsgc_
       basalmiddle.mod ihvcn.mod inav11.mod jsrnaf.mod ka.mod kcnq.mod kdpyr.mod kht.mod kif.mod 
       kir.mod kis.mod klt.mod kpkj.mod kpkj2.mod kpkjslow.mod kpksk.mod leak.mod multisite.mod 
       na.mod nacn.mod nacncoop.mod nap.mod napyr.mod pkjlk.mod rsg.mod vecevent.mod

     Testing with flags: -v --tb=short cnmodel/
    ============================= test session starts =============================
    platform win32 -- Python 2.7.14, pytest-3.3.2, py-1.5.2, pluggy-0.6.0 -- C:\Users\pbmanis\Anaconda2\python.exe
    cachedir: .cache
    rootdir: C:\Users\pbmanis\Desktop\Python\cnmodel, inifile:
    collecting ... collected 36 items

    cnmodel/an_model/tests/test_cache.py::test_cache PASSED                  [  2%]
    cnmodel/an_model/tests/test_cache.py::test_parallel PASSED               [  5%]
    cnmodel/cells/tests/test_cells.py::test_bushy PASSED                     [  8%]
    cnmodel/cells/tests/test_cells.py::test_bushy21 PASSED                   [ 11%]
    cnmodel/cells/tests/test_cells.py::test_bushy_mouse PASSED               [ 13%]
    cnmodel/cells/tests/test_cells.py::test_tstellate PASSED                 [ 16%]
    cnmodel/cells/tests/test_cells.py::test_tstellate_mouse PASSED           [ 19%]
    cnmodel/cells/tests/test_cells.py::test_tstellatet PASSED                [ 22%]
    cnmodel/cells/tests/test_cells.py::test_dstellate PASSED                 [ 25%]
    cnmodel/cells/tests/test_cells.py::test_dstellate_mouse PASSED           [ 27%]
    cnmodel/cells/tests/test_cells.py::test_octopus PASSED                   [ 30%]
    cnmodel/cells/tests/test_cells.py::test_pyramidal PASSED                 [ 33%]
    cnmodel/cells/tests/test_cells.py::test_tuberculoventral PASSED          [ 36%]
    cnmodel/cells/tests/test_cells.py::test_cartwheel PASSED                 [ 38%]
    cnmodel/cells/tests/test_cells.py::test_sgc_basal_middle PASSED          [ 41%]
    cnmodel/cells/tests/test_cells.py::test_sgc_apical PASSED                [ 44%]
    cnmodel/data/tests/test_db.py::test_db PASSED                            [ 47%]
    cnmodel/mechanisms/tests/test_mechanisms.py::test_max_open_probability PASSED [ 50%]
    cnmodel/synapses/tests/test_psd.py::test_sgc_bushy_psd PASSED            [ 52%]
    cnmodel/synapses/tests/test_psd.py::test_sgc_tstellate_psd PASSED        [ 55%]
    cnmodel/synapses/tests/test_psd.py::test_sgc_dstellate_psd PASSED        [ 58%]
    cnmodel/synapses/tests/test_psd.py::test_sgc_octopus_psd PASSED          [ 61%]
    cnmodel/synapses/tests/test_synapses.py::test_sgc_bushy PASSED           [ 63%]
    cnmodel/synapses/tests/test_synapses.py::test_sgc_tstellate PASSED       [ 66%]
    cnmodel/synapses/tests/test_synapses.py::test_sgc_tstellate2 PASSED      [ 69%]
    cnmodel/synapses/tests/test_synapses.py::test_sgc_dstellate PASSED       [ 72%]
    cnmodel/synapses/tests/test_synapses.py::test_dstellate_bushy PASSED     [ 75%]
    cnmodel/synapses/tests/test_synapses.py::test_dstellate_tstellate PASSED [ 77%]
    cnmodel/synapses/tests/test_synapses.py::test_dstellate_dstellate FAILED [ 80%]
    cnmodel/util/tests/test_expfitting.py::test_fit1 PASSED                  [ 83%]
    cnmodel/util/tests/test_expfitting.py::test_fit2 PASSED                  [ 86%]
    cnmodel/util/tests/test_matlab.py::test_matlab SKIPPED                   [ 88%]
    cnmodel/util/tests/test_sound.py::test_conversions PASSED                [ 91%]
    cnmodel/util/tests/test_sound.py::test_tonepip PASSED                    [ 94%]
    cnmodel/util/tests/test_sound.py::test_noisepip PASSED                   [ 97%]
    cnmodel/util/tests/test_stim.py::test_make_pulse PASSED                  [100%]

    ================================== FAILURES ===================================
    __________________________ test_dstellate_dstellate ___________________________
    cnmodel\synapses\tests\test_synapses.py:40: in test_dstellate_dstellate
        SynapseTester('dstellate', 'dstellate')
    cnmodel\synapses\tests\test_synapses.py:72: in __init__
        UserTester.__init__(self, "%s_%s" % (pre, post), pre, post)
    cnmodel\util\user_tester.py:33: in __init__
        self.assert_test_info(*args, **kwds)
    cnmodel\synapses\tests\test_synapses.py:108: in assert_test_info
        super(SynapseTester, self).assert_test_info(*args, **kwds)
    cnmodel\util\user_tester.py:127: in assert_test_info
        self.compare_results(result, expect)
    cnmodel\util\user_tester.py:60: in compare_results
        self.compare_results(info[k], expect[k])
    cnmodel\util\user_tester.py:63: in compare_results
        self.compare_results(info[i], expect[i])
    cnmodel\util\user_tester.py:79: in compare_results
        self.compare_results(info[k], expect[k])
    cnmodel\util\user_tester.py:71: in compare_results
        assert np.all(inans == enans)
    E   AssertionError
    ---------------------------- Captured stdout call -----------------------------
    << D-stellate: JSR Stellate Type I-II cell model created >>
    << D-stellate: JSR Stellate Type I-II cell model created >>
    Elapsed time for 1 Repetions: 0.355309
    =============== 1 failed, 34 passed, 1 skipped in 96.33 seconds ===============



Mac OSX (24 Feb 2018; Qt5 Branch)
---------------------

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




