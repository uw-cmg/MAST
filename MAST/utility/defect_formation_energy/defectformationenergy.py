import sys, os

from pymatgen.io.vaspio.vasp_output import Vasprun
from pymatgen.io.smartio import read_structure

from MAST.utility import PickleManager

class DefectFormationEnergy:
    def __init__(self, directory):
        self.directory = directory

        pm = PickleManager(self.directory + 'input_options.pickle')
        self.input_options = pm.load_variable()

    def calculate_defect_formation_energies(self):
        perf_dir = self.get_perfect_directory()[0]
        def_dir = self.get_defect_directories()

        defects = self.input_options.get_item('defects', 'defects')
        chempot = self.input_options.get_item('chemical_potentials')

        # Need to reconstitute the original system name
        system_name = self.input_options.get_item('mast', 'system_name').split('_')
        name = system_name[0]
        for string in system_name[1:-1]:
            name += ('_' + string)

        #print perf_dir
        print defects
        print chempot
        print name


        e_perf = self.get_total_energy(perf_dir)
        struct_perf = self.get_structure(perf_dir)
        perf_species = dict()
        for site in struct_perf:
            if (str(site.specie) not in perf_species):
                perf_species[str(site.specie)] = 1
            else:
                perf_species[str(site.specie)] += 1

        print 'Perfect composition: ', perf_species
        e_defects = dict() # dict of defect energies

        for label, defect in defects.items():
            def_dir = '%s/%s_%s_sp' % (self.directory, name, label)
            defect_info = dict()

            energy = self.get_total_energy(def_dir)
            structure = self.get_structure(def_dir)
            def_species = dict()
            for site in structure.sites:
                if (site.specie not in def_species):
                    def_species[site.specie] = 1
                else:
                    def_species[site.specie] += 1

            struct_diff = dict()
            for specie, number in def_species.items():
                try:
                    nperf = perf_species[str(specie)]
                except KeyError:
                    nperf = 0
                #print number - nperf
                struct_diff[str(specie)] = number - nperf

            base = energy - e_perf # E_defect - E_perf

            defect_info = dict()
            for conditions, potentials in chempot.items():
                print 'DFE for defect %s under %s conditions' % (label, conditions.upper())
                print '%-20s = %18.6f' % ('Energy defected', energy)
                print '%-20s = %18.6f' % ('Energy perfect', e_perf)
                e_def = base
                for specie, number in struct_diff.items():
                    mu = potentials[str(specie).lower()]
                    #print str(specie), mu, number
                    e_def -= (number * mu)
                defect_info[conditions] = e_def
                print 'DFE =', e_def

            e_defects[label] = defect_info
        return e_defects

    def get_defect_directories(self):
        return [x[0] for x in os.walk(self.directory) if ('defect_' in x[0] and '_sp' in x[0])]

    def get_perfect_directory(self):
        return [x[0] for x in os.walk(self.directory) if ('perfect' in x[0] and '_sp' in x[0])]

    def get_total_energy(self, directory):
        """Returns the total energy from a directory"""
        # Are we dealing with VASP?
        if ('vasprun.xml' in os.listdir(directory)):
        # Modified from the PyMatGen Vasprun.final_energy() function to return E_0
            return Vasprun('%s/vasprun.xml' % directory).ionic_steps[-1]["electronic_steps"][-1]["e_0_energy"]

    def get_structure(self, directory):
        """Returns the final structure from an optimization"""
        if ('vasprun.xml' in os.listdir(directory)):
            return Vasprun('%s/vasprun.xml' % directory).final_structure
