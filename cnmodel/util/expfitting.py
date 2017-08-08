#!/usr/bin/env python
# encoding: utf-8
"""
expfitting.py
Provide single or double exponential fits to data.
"""

import lmfit
import numpy as np
import scipy.optimize

class ExpFitting:
    """
    Parameters
    ----------
    nexp : int
        1 or 2 for single or double exponential fit
    initpars : dict
        dict of initial parameters. For example: {'dc': 0., 
        'a1': 1., 't1': 3, 'a2' : 0.5, 'delta': 3.}, where
        delta determines the ratio between the time constants.
    bounds : dict
        dictionary of bounds for each parameter, with a list of lower and upper values.
    
    """
    def __init__(self, nexp=1, initpars=None, bounds=None):
        self.fitpars = lmfit.Parameters()
        if nexp == 1:
            # (Name,  Value,  Vary,   Min,  Max,  Expr)
            self.fitpars.add_many(('dc', 0, True, -100., 0., None),
                                  ('a1', 1., True, -25., 25., None),
                                  ('t1', 10., True, 0.1, 50, None))
            self.efunc = self.exp1_err
        elif nexp == 2:
            self.fitpars.add_many(('dc', 0, True, -100., 0., None),
                                  ('a1', 1., True, 0., 25., None),
                                  ('t1', 10., True, 0.1, 50, None),
                                  ('a2', 1., True, 0., 25., None),
                                  ('delta', 3., True, 3., 100., None))
            if initpars is not None:
                assert len(initpars) == 5
                for k, v in initpars.iteritems():
                    self.fitpars[k].value = v
            if bounds is not None:
                assert len(bounds) == 5
                for k, v in bounds.iteritems():
                    self.fitpars[k].min = v[0]
                    self.fitpars[k].max = v[1]
                    
            self.efunc = self.exp2_err
        else:
            raise ValueError

    def fit(self, x, y, p, verbose=False):
        
        kws={'maxfev': 5000}
        mim = lmfit.minimize(self.efunc, p, method='least_squares', args=(x, y)) #, kws=kws)
        if verbose:
            lmfit.printfuncs.report_fit(mim.params)
        fitpars = mim.params
        return fitpars

    @staticmethod
    def exp1(x, dc, t1, a1):
        return dc + a1*np.exp(-x/t1)

    def exp1_err(self, p, x, y):
        return np.fabs(y-self.exp1(x, **dict([(k,p.value) for k,p in p.items()])))

    @staticmethod
    def exp2(x, dc, t1, a1, a2, delta):
        return dc + a1 * np.exp(-x/t1) + a2 * np.exp(-x/(t1*delta))

    def exp2_err(self, p, x, y):
        return np.fabs(y-self.exp2(x, **dict([(k,p.value) for k,p in p.items()])))


