#!/usr/bin/env python
__author__ = "Paul B. Manis"
__version__ = "0.2"

#import util
try:
    import faulthandler
    faulthandler.enable()
except ImportError:
    pass

import os
dirname = os.path.abspath(os.path.dirname(__file__))
libpath = os.path.join(dirname, '..')
import neuron
try:
    neuron.h.MultiSiteSynapse
except AttributeError:
    neuron.load_mechanisms(libpath)
# flag to allow unit tests to store / overwrite test results
AUDIT_TESTS = False



