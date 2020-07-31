Install of cnmodel 
------------------

cnmodel can be installed in several ways. 

The most "closed" and reproducible way is to use the
existing tools in python to create a virtual environment.

    1. Do a standard installation of NEURON

    2. Create virtual environment:

    (This scripts is in make_env.sh)

        ENVNAME="cnmodel_venv"
        python3 -m venv $ENVNAME
        source $ENVNAME/bin/activate
        pip install --upgrade pip  # be sure pip is up to date in the new env.
        pip install wheel  # seems to be missing (note singular)
        pip install cython
        # if requirements.txt is not present, you need to create it.:
        # However, we provide this in the files
        # pip install pipreqs
        # pipreqs
        #
        # Then:
        #
        pip install -r requirements.txt
        source $ENVNAME/bin/activate

        # build the mechanisms
        # this may equire a separate install of the standard NEURON package
        # with the same version as we have provided
        nrnivmodl cnmodel/mechanisms

        # this has to be done separately as the compilation depends on what happens above
        pip3 install -e git+https://github.com/pbmanis/cochlea-1.git@c2e8c9612481ebe397ba0e7762b8f2772500388d#egg=cochlea

        source $ENVNAME/bin/activate
        python setup.py develop

        # note that you will need to activate the environment once this script exists.
    
    3. in .zshrc, you can create a convenience alias to switch the environment and get into
       the directory.
   
ALternate:
    Install using anaconda python, building all by hand.
    