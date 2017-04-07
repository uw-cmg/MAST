#!/usr/bin/env python
######
# TTM 2017-04-07 convert old style XDATCAR with blank lines and Konfig lines
#                to Direct configuration
#####
import os
import sys

def general_instr():
    instr=list()
    instr.append("Add lattice parameters at the top.")
    instr.append("Top should then have formatting similar to the following:")
    instr.append("")
    instr.append("unknown system")
    instr.append("           1")
    instr.append("     7.878223    0.000048   -0.018810")
    instr.append("    -0.000087    7.871433    0.000077")
    instr.append("    -0.018798   -0.000063    7.878071")
    instr.append("   La   Ga   O ")
    instr.append("   8   8  23")
    instr.append("Direct configuration=     1")
    instr.append("")
    for myline in instr:
        print(myline)
    return

def replace_blanks(fname):
    with open(fname, 'r') as xfile:
        xlines = xfile.readlines()
    blankidxlist=list()
    for xidx in range(0, len(xlines)):
        xline = xlines[xidx]
        if len(xline.strip()) == 0:
            blankidxlist.append(xidx)
        elif "Konfig" in xline:
            blankidxlist.append(xidx)
    print(blankidxlist)
    for bidx in range(0, len(blankidxlist)):
        replace_idx = blankidxlist[bidx]
        confignum = bidx + 1 #start numbering at 1
        newline = "Direct configuration=" + ("%i" % confignum).rjust(6) + "\n"
        xlines[replace_idx] = newline
    fname_converted = fname + "_converted"
    with open(fname_converted, 'w') as xfile_converted:
        for xline in xlines:
            xfile_converted.write(xline)
    return fname_converted

def main(fname="XDATCAR"):
    print("Replacing blank lines.")
    fname_converted = replace_blanks(fname)
    print("Done in %s" % fname_converted)
    print("Done.")
    general_instr()
    return


if __name__=="__main__":
    fname="XDATCAR"
    if len(sys.argv) > 1:
        fname = sys.argv[1]
    main(fname)
