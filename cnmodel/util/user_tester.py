from __future__ import print_function
import os, sys, pickle, pprint
import numpy as np
import pyqtgraph as pg
from .. import AUDIT_TESTS

class UserTester(object):
    """
    Base class for testing when a human is required to verify the results.
    
    When a test is passed by the user, its output is saved and used as a basis
    for future tests. If future test results do not match the stored results,
    then the user is asked to decide whether to fail the test, or pass the
    test and store new results.
    
    Subclasses must reimplement run_test() to return a dictionary of results
    to store. Optionally, compare_results and audit_result may also be
    reimplemented to customize the testing behavior.
    
    By default, test results are stored in a 'test_data' directory relative
    to the file that defines the UserTester subclass in use.
    """
    data_dir = 'test_data'
    
    def __init__(self, key, *args, **kwds):
        """Initialize with a string *key* that provides a short, unique 
        description of this test. All other arguments are passed to run_test().
        
        *key* is used to determine the file name for storing test results. 
        """
        self.audit = AUDIT_TESTS
        self.key = key
        self.rtol = 1e-3
        self.args = args
        self.assert_test_info(*args, **kwds)
    
    def run_test(self, *args, **kwds):
        """
        Exceute the test. All arguments are taken from __init__.
        Return a picklable dictionary of test results.
        """
        raise NotImplementedError()

    def compare_results(self, key, info, expect):
        """
        Compare *result* of the current test against the previously stored 
        result *expect*. If *expect* is None, then no previous result was 
        stored.
        
        If *result* and *expect* do not match, then raise an exception.
        """
        # Check test structures are the same
        assert type(info) is type(expect)
        if hasattr(info, '__len__'):
            if len(info) != len(expect):
                print('\nComparing Lists, Model tested: %s, measure: %s' % (self.key, key))
                print('Lengths: expected: %d  received: %d' % (len(expect), len(info)))
                print('Array received: ', info)
                print('Array expected: ', expect)
            assert len(info) == len(expect)
            
        if isinstance(info, dict):
            for k in info:
                assert k in expect
            for k in expect:
                assert k in info
                self.compare_results(k, info[k], expect[k])
        elif isinstance(info, list):
            for i in range(len(info)):
                self.compare_results(key, info[i], expect[i])
        elif isinstance(info, np.ndarray):
            assert info.shape == expect.shape
            if len(info) == 0:
                return
            #assert info.dtype == expect.dtype
            if info.dtype.fields is None:
                intnan = -9223372036854775808  # happens when np.nan is cast to int
                inans = np.isnan(info) | (info == intnan)
                enans = np.isnan(expect) | (expect == intnan)
                assert np.all(inans == enans)
                mask = ~inans
                if not np.allclose(info[mask], expect[mask], rtol=self.rtol):
                    print( '\nComparing data array, shapes match: ', info.shape == expect.shape)
                    print( 'Model tested: %s, measure: %s' % (self.key, key))
                    print( 'args: ', dir(self.args[0]))
                    print( 'Array expected: ', expect[mask])
                    print('Array received: '  , info[mask])
                    print('args[0]: ', self.args[0])
                    #self.args[0].print_all_mechs()
                    try:
                        self.args[0].print_all_mechs()
                    except:
                        print('args[0] is string: ', self.args[0])
                assert np.allclose(info[mask], expect[mask], rtol=self.rtol)
            else:
                for k in info.dtype.fields.keys():
                    self.compare_results(k, info[k], expect[k])
        elif np.isscalar(info):
            if not np.allclose(info, expect, rtol=self.rtol):
                print('\nComparing Scalar data, model: %s, measure: %s' % (self.key, key))
                #print 'args: ', dir(self.args[0])
                print('Expected: ', expect, ',  received: ', info, '  relative tolerance: ', self.rtol)
                print('Ratio: received/expected: ', info/expect)
                print('Difference: received-expected: ', info - expect)
                if isinstance(self.args[0], str):
                    pass
                    # print ': ', str
                else:
                    self.args[0].print_all_mechs()
            assert np.allclose(info, expect, rtol=self.rtol)
        else:
            try:
                assert info == expect
            except AssertionError:
                raise
            except Exception:
                raise NotImplementedError("Cannot compare objects of type %s" % type(info))

    def audit_result(self, info, expect):
        """ Display results and ask the user to decide whether the test passed.
        Return True for pass, False for fail.
        
        If *expect* is None, then no previous test results were stored.
        """
        app = pg.mkQApp()
        print("\n=== New test results for %s: ===\n" % self.key)
        pprint.pprint(info)
        
        # we use DiffTreeWidget to display differences between large data structures, but
        # this is not present in mainline pyqtgraph yet
        if hasattr(pg, 'DiffTreeWidget'):
            win = pg.DiffTreeWidget()
        else:
            from cnmodel.util.difftreewidget import DiffTreeWidget
            win = DiffTreeWidget()
        
        win.resize(800, 800)
        win.setData(expect, info)
        win.show()
        print("Store new test results? [y/n]",)
        yn = input()
        win.hide()
        return yn.lower().startswith('y')
    
    def assert_test_info(self, *args, **kwds):
        """
        Test *cell* and raise exception if the results do not match prior
        data.
        """
        result = self.run_test(*args, **kwds)
        expect = self.load_test_result()
        try:
            assert expect is not None
            self.compare_results(None, result, expect)
        except:
            if not self.audit:
                if expect is None:
                    raise Exception("No prior test results for test '%s'. "
                                    "Run test.py --audit store new test data." % self.key)
                else:
                    raise
                
            store = self.audit_result(result, expect)
            if store:
                self.save_test_result(result)
            else:
                raise Exception("Rejected test results for '%s'" % self.key)
                    
    
    def result_file(self):
        """
        Return a file name to be used for storing / retrieving test results
        given *self.key*.
        """
        modfile = sys.modules[self.__class__.__module__].__file__
        path = os.path.dirname(modfile)
        return os.path.join(path, self.data_dir, self.key + '.pk')

    def load_test_result(self):
        """
        Load prior test results for *self.key*.
        If there are no prior results, return None.
        """
        fn = self.result_file()
        if os.path.isfile(fn):
            return pickle.load(open(fn, 'rb'),encoding='latin1')
        return None

    def save_test_result(self, result):
        """
        Store test results for *self.key*.
        Th e*result* argument must be picklable.
        """
        fn = self.result_file()
        dirname = os.path.dirname(fn)
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        pickle.dump(result, open(fn, 'wb'))
