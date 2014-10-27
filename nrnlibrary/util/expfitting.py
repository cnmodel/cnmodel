#!/usr/bin/env python
# encoding: utf-8
"""
expfitting.py

Created by Paul Manis on 2014-10-08.

"""

import sys
import os
import unittest
import lmfit
import numpy as np
#           (Name,  Value,  Vary,   Min,  Max,  Expr)
# p.add_many(('amp1',    10,  True, None, None,  None),
#            ('cen1',   1.2,  True,  0.5,  2.0,  None),
#            ('wid1',   0.8,  True,  0.1, None,  None),
#            ('amp2',   7.5,  True, None, None,  None),
#            ('cen2',   1.9,  True,  1.0,  3.0,  None),
#            ('wid2',  None, False, None, None, '2*wid1/3'))
class ExpFitting:
    def __init__(self, nexp=1, initpars=None, bounds=None):
        self.fitpars = lmfit.Parameters()
        if nexp == 1:
            # (Name,  Value,  Vary,   Min,  Max,  Expr)
            self.fitpars.add_many(('dc', 0, True, -100., 0., None),
                                  ('a1', 1., True, -25., 25., None),
                                  ('t1', 10., True, 0.1, 50, None))
            self.efunc = self.exp1
        elif nexp == 2:
            self.fitpars.add_many(('dc', 0, True, -100., 0., None),
                                  ('a1', 1., True, 0., 25., None),
                                  ('t1', 10., True, 0.1, 50, None),
                                  ('a2', 1., True, 0., 25., None),
                                  ('delta', 3., True, 3., 100., None))
            if initpars is not None and len(initpars) == 5:
                for k, v in initpars.iteritems():
                    self.fitpars[k].value = v
            if bounds is not None and len(bounds) == 5:
                for k, v in bounds.iteritems():
                    self.fitpars[k].min = v[0]
                    self.fitpars[k].max = v[1]
                    
            self.efunc = self.exp2
        else:
            raise ValueError

    def fit(self, x, y, p, verbose=False):
        
        mim = lmfit.minimize(self.efunc, p, args=(x, y))
        mim.leastsq()
        if verbose:
            lmfit.printfuncs.report_fit(mim.params)
        fitpars = mim.params
        return fitpars

    def exp1(self, p, x, y):
        yd = p['dc'].value + p['a1'].value*np.exp(-x/p['t1'].value)
        return np.fabs(y-yd)

    def exp2(self, p, x, y):
        yd = p['dc'].value + p['a1'].value*np.exp(-x/p['t1'].value) + p['a2'].value*np.exp(-x/(p['t1'].value*p['delta'].value))
        return np.fabs(y-yd)
        

class TestsExpFitting(unittest.TestCase):
    
    def setUp(self):
        pass
        
    def test_fit(self):
        self.fit = ExpFitting(nexp=1)
        self.x = np.linspace(0., 50, 500)
        self.p = [-50., 4., 5.]
        self.y = self.p[0] + self.p[1]*np.exp(-self.x/self.p[2])
        res = self.fit.fit(self.x, self.y, self.fit.fitpars)
        pr = [float(res[k].value) for k in res.keys()]
        print '\noriginal: ', self.p
        print 'fit res:  ', pr
        for i, v in enumerate(self.p):
            self.assertAlmostEqual(v, pr[i], delta=1e-6)

    def test_fit2(self):
        self.fit = ExpFitting(nexp=2)
        self.x = np.linspace(0., 50, 500)
        self.p = [-50., 4., 5., 1., 4.5]   # last term is ratio of the two time constants (t2 = delta*t1)
        self.y = self.p[0] + self.p[1]*np.exp(-self.x/self.p[2]) + self.p[3]*np.exp(-self.x/(self.p[2]*self.p[4]))
        res = self.fit.fit(self.x, self.y, self.fit.fitpars)
        pr = [float(res[k].value) for k in res.keys()]
        print '\noriginal: ', self.p
        print 'fit res:  ', pr
        # we can only do this approximately for 2 exp fits
        for i, v in enumerate(self.p): # test each one individually
            self.assertAlmostEqual(pr[i]/v, 1.0 , delta=1e-4)
        

if __name__ == '__main__':
    #unittest.main()
    print 'entering tests'
    suite = unittest.TestLoader().loadTestsFromTestCase(TestsExpFitting)
    print 'suite: ', suite
    unittest.TextTestRunner(verbosity=2).run(suite)