"""
Run unit tests for nrnlibrary
"""

import os, sys
import pytest

# Make sure we look for nrnlibrary here first.
path = os.path.dirname(__file__)
sys.path.insert(0, path)


# Allow user to audit tests with --audit flag
import nrnlibrary
if '--audit' in sys.argv:
    sys.argv.remove('--audit')
    sys.argv.append('-s') # needed for cli-based user interaction
    nrnlibrary.AUDIT_TESTS = True

# generate test flags
flags = sys.argv[1:]
tb = [flag for flag in flags if flag.startswith('--tb')]
if len(tb) == 0:
    flags.append('--tb=short')

add_path = True
for flag in flags:
    if os.path.isdir(flag) or os.path.isfile(flag):
        add_path = False
        break
if add_path:
    flags.append('nrnlibrary/')

# Start tests.
print "Testing with flags: %s" % " ".join(flags)
pytest.main(flags)
