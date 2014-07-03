#-----------------------------------------------------------------------------
# Standard Hodkgin-Huxley model, small cell


def hh(debug=False, message=None):
    """
    Standard Hodgkin-Huxley mechanisms from NEURON
    """
    soma = h.Section() # one compartment of about 29000 um2
    v_potassium = -80       # potassium reversal potential
    v_sodium = 50           # sodium reversal potential
    c_m = 1.0
    scalefactor = 1.0 # This determines the relative size of the cell
    rinsf = 1.0           # input resistance adjustment (also current...)
    totcap = 20.0 # scalefactor * 1.0 # cap in pF for cell
    effcap = totcap # sometimes we change capacitance - that's effcap
    somaarea = totcap * 1E-6 / c_m # pf -> uF, cm = 1uf/cm^2 nominal
    lstd = 1E4 * sqrt(somaarea / 3.14159) # convert from cm to um

    soma.nseg = 1
    soma.diam = lstd
    soma.L = lstd

    seg = soma
    seg.insert('hh')
    seg.insert('pas')
    if not runQuiet:
        if message is None:
            print "<< Standard HH model created >>"
        else:
            print message
    return(soma, (None, None, None))
