import lmfit
import numpy as np


class FitModel(lmfit.Model):
    """ Simple extension of lmfit.Model that allows one-line fitting.
    
    Example uses:
        
        # single exponential fit
        fit = expfitting.Exp1.fit(data, 
                x=time_vals,
                xoffset=(0, 'fixed'),
                yoffset=(yoff_guess, -120, 0),
                amp=(amp_guess, 0, 50),
                tau=(tau_guess, 0.1, 50))
        
        # plot the fit
        fit_curve = fit.eval()
        plot(time_vals, fit_curve)
        
    
        # double exponential fit with tau ratio constraint
        # note that 'tau_ratio' does not appear in the exp2 model; 
        # we can define new parameters here.
        fit = expfitting.Exp1.fit(data, 
                x=time_vals,
                xoffset=(0, 'fixed'),
                yoffset=(yoff_guess, -120, 0),
                amp1=(amp_guess, 0, 50),
                tau1=(tau_guess, 0.1, 50),
                amp2=(-0.5, -50, 0),
                tau_ratio=(10, 3, 50),
                tau2='tau1 * tau_ratio'
                )
        
    """
    def fit(self, data, **params):
        """ Return a fit of data to this model.
        
        Parameters
        ----------
        data : array
            dependent data to fit
            
        Extra keyword arguments are passed to make_params() if they are model
        parameter names, or passed directly to Model.fit() for independent
        variable names.
        """
        fit_params = {}
        model_params = {}
        for k,v in params.items():
            if k in self.independent_vars or k in ['weights', 'method', 'scale_covar', 'iter_cb']:
                fit_params[k] = v
            else:
                model_params[k] = v
        p = self.make_params(**model_params)
        return lmfit.Model.fit(self, data, params=p, **fit_params)
        
    def make_params(self, **params):
        """
        Make parameters used for fitting with this model.
        
        Keyword arguments are used to generate parameters for the fit. Each 
        parameter may be specified by the following formats:
        
        param=value :
            The initial value of the parameter
        param=(value, 'fixed') :
            Fixed value for the parameter
        param=(value, min, max) :
            Initial value and min, max values, which may be float or None
        param='expression' :
            Expression used to compute parameter value. See:
            http://lmfit.github.io/lmfit-py/constraints.html#constraints-chapter
        """
        p = lmfit.Parameters()
        for k in self.param_names:
            p.add(k)
        
        for param,val in params.items():
            if param not in p:
                p.add(param)
                
            if isinstance(val, str):
                p[param].expr = val
            elif np.isscalar(val):
                p[param].value = val
            elif isinstance(val, tuple):
                if len(val) == 2:
                    assert val[1] == 'fixed'
                    p[param].value = val[0]
                    p[param].vary = False
                elif len(val) == 3:
                    p[param].value = val[0]
                    p[param].min = val[1]
                    p[param].max = val[2]
                else:
                    raise TypeError("Tuple parameter specifications must be (val, 'fixed')"
                                    " or (val, min, max).")
            else:
                raise TypeError("Invalid parameter specification: %r" % val)
            
        # set initial values for parameters with mathematical constraints
        # this is to allow fit.eval(**fit.init_params)
        global_ns = np.__dict__
        for param,val in params.items():
            if isinstance(val, str):
                p[param].value = eval(val, global_ns, p.valuesdict())
        return p
        
        
def exp1(x, xoffset, yoffset, tau, amp):
    return yoffset + amp * np.exp(-(x - xoffset)/tau)

class Exp1(FitModel):
    """ Single exponential fitting model.
    
    Parameters are xoffset, yoffset, amp, and tau.
    """
    def __init__(self):
        FitModel.__init__(self, exp1, independent_vars=['x'])


def exp2(x, xoffset, yoffset, tau1, amp1, tau2, amp2):
    xoff = x - xoffset
    return yoffset + amp1 * np.exp(-xoff/tau1) + amp2 * np.exp(-xoff/tau2)

class Exp2(FitModel):
    """ Double exponential fitting model.
    
    Parameters are xoffset, yoffset, amp1, tau1, amp2, and tau2.
    """
    def __init__(self):
        FitModel.__init__(self, exp2, independent_vars=['x'])


