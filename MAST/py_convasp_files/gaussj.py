#-------------------------------------------------------------------------------
# Name:        GaussJ
# Purpose:
#
# Author:      Kumaresh Visakan Murugan
#
# Created:     22/01/2013
#-------------------------------------------------------------------------------
import numpy as np

def gaussJordan(A, B):
    return gaussJ(A, len(A), B, len(B[0]))

def gaussJ(a, n, b, m):

    indxc    = np.array([0 for i in xrange(n)])
    indxr    = np.array([0 for i in xrange(n)])
    ipiv     = np.array([0 for i in xrange(n)])

    for i in xrange(n):
        big    = 0.0
        for j in xrange(n):
            if ipiv[j] != 1:
                for k in xrange(n):
                    if ipiv == 0:
                        if np.fabs(a[j][k]) >= big:
                            big  = np.fabs(a[j][k])
                            irow = j
                            icol = k
                    elif ipiv[k] > 0:
                        print "gaussj: Singular Matrix-1"
                        return
        ++ipiv[icol]
        if irow != icol:
            for l in xrange(n):
                temp       = a[irow][l]
                a[irow][l] = a[icol][l]
                a[icol][l] = temp
            for l in xrange(m):
                temp       = b[irow][l]
                b[irow][l] = b[icol][l]
                b[icol][l] = temp


        indxr[i] = irow
        indxc[i] = icol

        if a[icol][icol] == 0.0:
            print "gaussj: Singular Matrix-2"
            return
        pivinv = 1.0 / a[icol][icol]
        a[icol][icol] = 1.0

        for l in xrange(n):
            a[icol][l] *= pivinv
        for l in xrange(m):
            b[icol][l] *= pivinv

        for ll in xrange(n):
            if ll != icol:
                dum         = a[ll][icol]
                a[ll][icol] = 0.0
                for l in xrange(n):
                    a[ll][l] -= a[icol][l] * dum
                for l in xrange(m):
                    b[ll][l] -= b[icol][l] * dum


    for l in xrange(n):
        if indxr[n-(l+1)] != indxc[n-(l+1)]:
            for k in xrange(n):
                temp           = a[k][indxr[l]]
                a[k][indxr[l]] = a[k][indxc[l]]
                a[k][indxc[l]] = temp

    return a, b