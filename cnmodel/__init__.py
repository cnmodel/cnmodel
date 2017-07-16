__author__ = "Paul B. Manis and Luke Campagnola"
__version__ = "0.3"

try:
    import faulthandler
    faulthandler.enable()
except ImportError:
    pass

import logging
logging.basicConfig(level=logging.INFO, format="[%(process)s] %(message)s")
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



