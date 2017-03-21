#!/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2017-03-21
##############################################################
import optparse
import pymatgen
import dirutil
import numpy as np
from MAST.utility import MASTFile
def main():
    """Takes a Diffusivity file and some inputs and returns host-value
        adjusted values. See Wu et al., Scientific Data, 2016.
    """
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input', dest='input', default='Diffusivity.txt',
                      help="Diffusivity text file.")
    parser.add_option('-q', '--hostqadjust', dest='hostqadjust', default='0.0',
                      help="Additive Q adjustment for host system.")
    parser.add_option('-d', '--hostdadjust', dest='hostdadjust', default='1.0',
                      help="Multiplicative D0 adjustment for host system.")
    parser.add_option('-m', '--tmelt', dest='tmelt', default='1400',
                        help="Melting temperature (K).")
    
    (myopt, myarg) = parser.parse_args()
    print myopt
    print myarg

    print "---------------------------------------------------"
    print "Welcome to the MAterials Simulation Toolkit (MAST)"
    print "    Pymatgen symmetry-finder wrapper tool"
    print "Version: " + dirutil.get_version()
    print "Installed in: " + dirutil.get_mast_install_path()
    print "---------------------------------------------------"
   
    myfile = MASTFile(myopt.input)
    mylines = list(myfile.data)
    mylines.pop(0) # remove header line of 1000/T in K and ln(D)
    parsedlist=list()
    for myline in mylines:
        [thousand_over_T, my_D] = myline.strip().split()
        thousand_over_T = float(thousand_over_T)
        my_D = float(my_D)
        one_over_T = thousand_over_T / 1000.0
        ln_D = np.log(my_D)
        parsedlist.append([one_over_T, ln_D, thousand_over_T, my_D])
    tmelt = float(myopt.tmelt)
    tmax = tmelt * 1.1 #pad a little bit, since temperature spacing is inexact
    tmin = 0.4 * tmelt #pad a little bit from 0.5
    shortlist=list()
    for parseditem in parsedlist:
        one_over_T = parseditem[0]
        if one_over_T == 0.0:
            my_T = 1e20 #approximate infinity with unreasonably high K temp
        else:
            my_T = 1.0/one_over_T
        if (my_T <= tmax) and (my_T >= tmin):
            shortlist.append(parseditem)
    print("1/T and ln(D), and 1000/T and D within 0.5*Tmelt and Tmelt:")
    for shortitem in shortlist:
        print("%3.5E   %3.3E   %3.5E    %3.3E" % (shortitem[0],shortitem[1],
                    shortitem[2],shortitem[3]))
    shortarr = np.array(shortlist)
    xdata = shortarr[:,0]
    ydata = shortarr[:,1]
    [slope, intercept] = np.polyfit(xdata, ydata, 1)
    kboltz = 0.00008617 # eV/K
    print "Slope %s and intercept %s" % (slope, intercept)
    qunadjusted = -1 * slope * kboltz
    dunadjusted = np.exp(intercept)
    print("Undadjusted Q: %3.3f eV" % qunadjusted)
    print("Undadjusted D0: %3.3f cm^2/sec" % dunadjusted)
    hostqadjust = float(myopt.hostqadjust)
    hostdadjust = float(myopt.hostdadjust)
    qadjusted = qunadjusted + hostqadjust
    dadjusted = dunadjusted * hostdadjust
    print("Adjusted Q: %3.3f eV" % qadjusted)
    print("Adjusted D0: %3.3f cm^2/sec" % dadjusted)
    return [qadjusted, dadjusted]

if __name__ == '__main__':
    #Use example from command prompt:
    #python diffcoeff_adjustment_util.py -i <path_to_mast>/MAST/test/diffusion_coefficient_test/files/diffcoeff_utility_bcc_9freq_WAg/Diffusivity.txt -m 3695 -q 1.3 -d 800
    main()
