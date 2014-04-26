##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Glen Jenness
# Last updated: 2013-07-01
##############################################################

import sys, math
import numpy as np


class PotentialAlignment:
#    def __init__(self):
#        self.outcar_perf = outcar_perf
#        self.outcar_def = outcar_def

    def read_outcar(self, outcar):
        """Reads a VASP OUTCAR file, and extracts out some useful information.
            Args:
                outcar:  OUTCAR file name
            Returns:
                elst_pot:  average electrostatic potential for each atomic center
                atom_list:  list of each atom type
        """
        outfile = open(outcar, 'r')

        atom_types = list()
        atom_numbers = list()
        atom_list = list()

        eof = False
        while (not eof):
            line = outfile.readline()
            if (not line):
                eof = True
            elif 'average (electrostatic) potential at core' in line:
                raw_data = list()
                outfile.readline() # Read in Junk line
                outfile.readline() # Read in Junk line

                pot_line = True
                while (pot_line):
                    pot_line = outfile.readline().strip()
                    if (not line):
                        pot_line = False
                    else:
                        raw_data.append(pot_line)
            elif 'ions per type' in line:
                atom_numbers = map(int, line.split('=')[1].split())
            elif 'VRHFIN' in line:
                atom_types.append(line.split()[1][1:-1])
            elif ' E-fermi :' in line:
                efermi = float(line.split()[2])

        outfile.close()
        elst_pot = list()
#        print raw_data
        for data in raw_data:
            data = data.split('-')
            for datum in data[1:-1]:
                elst_pot.append(-float(datum.split()[0]))
            if data[0]:
                elst_pot.append(-float(data[-1]))

        for i in range(len(atom_numbers)):
            for j in range(atom_numbers[i]):
                atom_list.append(atom_types[i])

        return (elst_pot, atom_list, efermi)

    def get_potential_alignment(self, perfect_info, defect_info):
        """Gets the shift.  Assumes all atoms in the perfect cell are symmetry
            equivalent, i.e. for Sc2O3, by symmetry we have two different Sc
            atoms (8b and 24d Wycoff positions) with two different potentials.

            In reality the best way to do this is to match up the symmetry
            positions of the perfect with that in the defected and do the 
            alignment that way.  However, a quick and dirty way is to average
            over each symmetry position, which allows us to account for the 
            different symmetry environments in an effective way without having
            to dig into the underlying unit cell symmetry.
        """

        perf_elst_pot, perf_atom_list, perf_efermi = perfect_info
        def_elst_pot, def_atom_list, def_efermi = defect_info

        unique_perf = list(set(perf_atom_list))
        ave_pot = list()

        for unique in unique_perf:
            nunique = 0
            unique_sum = 0.0
            for i in range(len(perf_atom_list)):
                if (perf_atom_list[i] == unique):
                    unique_sum += perf_elst_pot[i]
                    nunique += 1

            ave_pot.append(unique_sum / nunique)

#        print ave_pot
        shift = 0.0

        for i in range(len(ave_pot)):
            perf_pot = ave_pot[i]
            atom = unique_perf[i]

            for j in range(len(def_elst_pot)):
                if (atom == def_atom_list[j]):
                    shift += (def_elst_pot[j] - perf_pot)

#        print perf_efermi
        return shift / len(def_elst_pot)
