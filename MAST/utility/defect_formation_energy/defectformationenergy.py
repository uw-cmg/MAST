import sys, os

from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.io.smartio import read_structure

from MAST.utility import PickleManager
from MAST.utility.defect_formation_energy.potential_alignment import PotentialAlignment

class DefectFormationEnergy:
    """Class for calculating the defect formation energy for a completed MAST
        run.
    """

    def __init__(self, directory):
        self.directory = directory

        pm = PickleManager(self.directory + '/input_options.pickle')
        self.input_options = pm.load_variable()
        self.recipe_plan = PickleManager(self.directory + '/mast.pickle').load_variable()

        self.final_ingredients = list()
        for ingredient in self.recipe_plan:
            if (ingredient.children is None):
                self.final_ingredients.append(ingredient.name)

    def calculate_defect_formation_energies(self):
        perf_dir = [ingredient for ingredient in self.final_ingredients if ('perfect' in ingredient)][0]
        def_dir = [ingredient for ingredient in self.final_ingredients if 'perfect' not in ingredient] 

        defects = self.input_options.get_item('defects', 'defects')
        chempot = self.input_options.get_item('chemical_potentials')

        # Need to reconstitute the original system name
        system_name = self.input_options.get_item('mast', 'system_name').split('_')
        name = system_name[0]
        for string in system_name[1:-1]:
            name += ('_' + string)

        e_perf = self.get_total_energy(perf_dir)
        efermi = self.get_fermi_energy(perf_dir)
        struct_perf = self.get_structure(perf_dir)

        perf_species = dict()
        for site in struct_perf:
            if (str(site.specie) not in perf_species):
                perf_species[str(site.specie)] = 1
            else:
                perf_species[str(site.specie)] += 1

        #print 'Perfect composition: ', perf_species
        e_defects = dict() # dict of defect energies

        # First loop through the conditions for the defects
        for conditions, potentials in chempot.items():
            print 'Conditions:  %s' % conditions.upper()
            defect_info = dict()

            # Loop through each defect
            for label, defect in defects.items():
                print 'Calculating DFEs for defect %s' % (label)

                # List of directories pertaining to a certain defect
                ddir = [ddir for ddir in def_dir if (label.split('defect_')[1] == ddir.split('defect_')[1].split('_q=')[0])]

                # Now loop through the charges for each defect
                energy_list = list()
                for charged in ddir:
                    charge = charged.split('_q=')[-1].split('_')[0]
                    sign = charge[0]

                    if (sign == 'p'):
                        charge = int(charge[1:])
                    elif (sign == 'n'):
                        charge = -1 * int(charge[1:])

                    energy = self.get_total_energy(charged)
                    structure = self.get_structure(charged)

                    # Find out how many atoms of each type are in the defects
                    def_species = dict()
                    for site in structure.sites:
                        if (site.specie not in def_species):
                            def_species[site.specie] = 1
                        else:
                            def_species[site.specie] += 1

                    # Find the differences in the number of each atom type
                    # between the perfect and the defect
                    struct_diff = dict()
                    for specie, number in def_species.items():
                        try:
                            nperf = perf_species[str(specie)]
                        except KeyError:
                            nperf = 0
                        struct_diff[str(specie)] = number - nperf

                    # Get the potential alignment correction
                    alignment = self.get_potential_alignment(perf_dir, charged)

                    # Calculate the base DFE energy
                    e_def = energy - e_perf # E_defect - E_perf
                    for specie, number in struct_diff.items():
                        mu = potentials[str(specie).lower()]
                        #print str(specie), mu, number
                        e_def -= (number * mu)
                    e_def += charge * (efermi + alignment) # Add in the shift here!
                    energy_list.append([charge, e_def])

                energy_list = sorted(energy_list, key=lambda energy: energy[0])
                defect_info[label] = energy_list

            e_defects[conditions] = defect_info

        return e_defects

    def get_total_energy(self, directory):
        """Returns the total energy from a directory"""

        abspath = '%s/%s/' % (self.directory, directory)

        if ('vasprun.xml' in os.listdir(abspath)):
        # Modified from the PyMatGen Vasprun.final_energy() function to return E_0
            return Vasprun('%s/vasprun.xml' % abspath).ionic_steps[-1]["electronic_steps"][-1]["e_0_energy"]

    def get_fermi_energy(self, directory):
        """Returns the Fermi energy from a directory"""
        abspath = '%s/%s/' % (self.directory, directory)

        if ('vasprun.xml' in os.listdir(abspath)):
            return Vasprun('%s/vasprun.xml' % abspath).efermi

    def get_structure(self, directory):
        """Returns the final structure from an optimization"""

        abspath = '%s/%s/' % (self.directory, directory)

        if ('vasprun.xml' in os.listdir(abspath)):
            return Vasprun('%s/vasprun.xml' % abspath).final_structure

    def get_potential_alignment(self, perf_dir, def_dir):
        """Returns the potential alignment correction used in charge defects"""

        abs_path_perf = '%s/%s/' % (self.directory, perf_dir)
        abs_path_def = '%s/%s/' % (self.directory, def_dir)

        pa = PotentialAlignment()

        if ('OUTCAR' in os.listdir(abs_path_perf)):
            perfect_info = pa.read_outcar('%s/%s' % (abs_path_perf, 'OUTCAR'))
            defect_info = pa.read_outcar('%s/%s' % (abs_path_def, 'OUTCAR'))

            return pa.get_potential_alignment(perfect_info, defect_info)

    def print_table(self):
        e_defects = self.calculate_defect_formation_energies()

        for conditions, defects in e_defects.items():
            print '\n\nDefect formation energies for %s conditions.\n' % conditions.upper()
            print '%-20s | %10s' % ('Defect', 'DFE')
            print '---------------------------------'
            for defect, energies in defects.items():
                for energy in energies:
                    print '%-14s(q=%2i) | %8.4f' % (defect, energy[0], energy[1])
                print str() # Add a blank line here
            print '---------------------------------'

