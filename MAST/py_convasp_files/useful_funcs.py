#-------------------------------------------------------------------------------
# Name:        UsefulFuncs
# Purpose:
#
# Author:      Kumaresh Visakan
#
# Created:     23/01/2013
#-------------------------------------------------------------------------------

import numpy as np

def VVdiff(a, b):
    size = min(len(a), len(b))
    c    = np.array([0.0 for i in xrange(size)])
    for i in xrange(size):
        c[i] = a[i] = b[i]

    return c

def VVsum(a, b):
    size = min(len(a), len(b))
    c    = np.array([0.0 for i in xrange(size)])
    for i in xrange(size):
        c[i] = a[i] + b[i]
    return c

def SVprod(s, b):
    return s * b

def vecD2C(lat, vd):
    vc = np.array([0.0, 0.0, 0.0])
    for ic in xrange(3):
        for jc in xrange(3):
            vc[ic] = vc[ic] + vd[jc] * lat[jc][ic]
    return vc