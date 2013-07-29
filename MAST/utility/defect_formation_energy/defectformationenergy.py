import sys, os

from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.io.smartio import read_structure

from MAST.utility import PickleManager
from MAST.utility.defect_formation_energy.potential_alignment import PotentialAlignment
from MAST.utility import Metadata

class DefectFormationEnergy:
    """Class for calculating the defect formation energy for a completed MAST
        run.
    """

    def __init__(self, directory=None, plot_threshold=0.01):
        self.directory = directory
        self.plot_threshold = plot_threshold

        pm = PickleManager(self.directory + '/input_options.pickle')
        self.input_options = pm.load_variable()
        self.recipe_plan = PickleManager(self.directory + '/mast.pickle').load_variable()

        self.final_ingredients = list()
        for ingredient in self.recipe_plan:
            if (ingredient.children is None):
                self.final_ingredients.append(ingredient.name)
        #print self.final_ingredients

    def calculate_defect_formation_energies(self):
        perf_dir = [ingredient for ingredient in self.final_ingredients if ('perfect' in ingredient)][0]
        def_dir = [ingredient for ingredient in self.final_ingredients if 'perfect' not in ingredient] 
        for whatever in sorted(def_dir):
            print whatever

        defects = self.input_options.get_item('defects', 'defects')
        chempot = self.input_options.get_item('chemical_potentials')

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
            for ddir in sorted(def_dir):
                def_meta = Metadata(metafile='%s/%s/metadata.txt' % (self.directory, ddir))
                label = def_meta.read_data('defect_label')
                charge = int(def_meta.read_data('charge'))
                energy = float(def_meta.read_data('energy'))
                structure = self.get_structure(ddir)

                print 'Calculating DFEs for defect %s with charge %3i.' % (label, charge)

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
                alignment = self.get_potential_alignment(perf_dir, ddir)

                # Calculate the base DFE energy
                e_def = energy - e_perf # E_defect - E_perf
                for specie, number in struct_diff.items():
                    mu = potentials[str(specie)]
                    #print str(specie), mu, number
                    e_def -= (number * mu)
                e_def += charge * (efermi + alignment) # Add in the shift here!
                print 'DFE = %f' % e_def

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

