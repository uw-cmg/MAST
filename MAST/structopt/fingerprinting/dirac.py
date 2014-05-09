import math

def dirac(x,a=0,sig=None):
    """Simple calculation of the dirac function
    Inputs:
        x = position for function evaluation
        a = dirac move
        sig = dirac spread
    Outputs:
        out = value of dirac function evaluated at x
    """
    if sig==None:
        if x != a:
            out=0
        else:
            out=1
    else:
        out=1/(sig*math.pi**0.5)*math.exp(-(x-a)**2/sig**2)
    return out

