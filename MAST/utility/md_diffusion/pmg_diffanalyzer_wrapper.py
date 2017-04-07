#!/usr/bin/env python
############
# Tam Mayeshiba 2017-03-16, 2017-04-07
# This file and pmg_diffanalyzer_input.txt just wrap around
# pymatgen's diffusion analyzer
# using its from_structures creator.
# See pymatgen's documentation at 
#      http://pymatgen.org/_modules/pymatgen/analysis/diffusion_analyzer.html
# pymatgen.analysis.diffusion_analyzer wrapper
############
import os
import sys
from pymatgen.io.vasp import Xdatcar
from pymatgen.analysis.diffusion_analyzer import DiffusionAnalyzer as DA
import numpy as np

def parse_input(ifilename=""):
    idict=dict()
    with open(ifilename,'r') as ifile:
        ilines = ifile.readlines()
    for iline in ilines:
        sline = iline.split("#")[0] #remove comments
        sline = sline.strip()
        if len(sline) == 0:
            continue
        ssplit = sline.split("=")
        skey = ssplit[0].strip()
        sval = ssplit[1].strip()
        if (sval.isalpha()) or (skey in ['filename','specie']):
            pass
        else:
            sval = float(sval)
        idict[skey] = sval
    for key, value in idict.items():
        print("%s: %s" % (key, value))
    return idict


def main(filename="XDATCAR", specie=None,
            temperature=298,
            time_step = 1,
            step_skip = 1000,
            smoothed = None,
            min_obs = 30,
            avg_nsteps = 1000,
            *args,**kwargs):
    #parse XDATCAR
    myxdatcar = Xdatcar(filename)
    #get structure list
    str_list = list(myxdatcar.structures)
    if specie == None:
        raise ValueError("Specie cannot be None")
    print "%i structures" % len(str_list)
    #run DA
    myDA = DA.from_structures(str_list, 
                        specie, temperature, time_step, step_skip,
                        smoothed=smoothed,
                        min_obs=min_obs,
                        avg_nsteps=avg_nsteps,
                        initial_disp=None,
                        initial_structure=None)
    myDA.export_msdt("output_msd.csv")
    summary = myDA.get_summary_dict()
    print(summary)
    print("Diffusion coefficient (cm^2/sec):%s" % summary['D'])
    with open("output_d_summary.txt",'w') as ofile:
        for key, value in summary.items():
            ofile.write("%s: %s\n" % (key, value))
    return summary

if __name__ == "__main__":
    ifilename = "pmg_diffanalyzer_input.txt"
    if len(sys.argv) > 1:
        ifilename = sys.argv[1]
    kwargs = parse_input(ifilename)
    main(**kwargs)
    sys.exit()
