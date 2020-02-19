Installation
============

Requirements
------------

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
        $ git checkout python3 #(if necessary)

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


Windows Notes
-------------

1. For more detailed information on setup in a Windows environment for Python 2.7, see the file Windows_setup.md. Thanks to Laurel Carney for prompting the generation of this set of instructions, and for identifying issues on Windows. A similar setup should work for Python 3.6+.

2. Manually compile the mex files for the Zilany et al model. In Matlab, go to the an_model/models folder, and use mexANmodel.m to compile the files. Then, add the an_model/model folder to the Matlab path, so that it can find the files when needed.

3. Under Windows, it may be best to use the standard Windows command terminal rather than the "bash" terminal provided by NEURON, at least to run the Python scripts.


MATLAB (R)
----------

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
