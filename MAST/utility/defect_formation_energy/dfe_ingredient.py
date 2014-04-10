import sys, os

import pymatgen as pmg
from pymatgen.io.vaspio.vasp_output import Vasprun, Outcar
from pymatgen.io.smartio import read_structure
from MAST.utility import MASTError
from MAST.utility import MASTFile
from MAST.utility.defect_formation_energy.defectformationenergy import DefectFormationEnergy
#from MAST.utility import PickleManager
from MAST.utility.defect_formation_energy.potential_alignment import PotentialAlignment
from MAST.utility.defect_formation_energy.gapplot import GapPlot
from MAST.utility import Metadata
from MAST.parsers import InputParser

class DefectFormationEnergyIngredient(DefectFormationEnergy):
    """Class for calculating the defect formation energy for a completed MAST
        run, as an ingredient whose input file contains
        dfe_label=perfect_label defected_label
        bandgap_lda_or_gga=<float>
        bandgap_hse_or_expt=<float>
        Attributes:
            self.ingdir <str>: Ingredient directory (should be the
                                current directory)
            self.recdir <str>: Recipe directory (one level up from self.ingdir)
            self.plot_threshold <float>: Plotting threshold
            self.input_options <InputOptions>: Input options of recipe
            self.dirs <dict of str>: Dictionary of perfect and defected
                                     directories, as:
                self.dirs['label1']['perfect'] = <str, directory>
                self.dirs['label1']['defected'] = <str, directory>
            self.bandgap_lda_or_gga <float>: LDA or GGA bandgap (usually
                                             underestimated)
            self.bandgap_hse_or_expt <float>: HSE or experimental bandgap
            self.e_defects <dict>: Defect formation energy dictionary
                self.e_defects[conditions][label] = <float>
    """

    def __init__(self, inputtxt="", plot_threshold=0.01):
        """Initialize defect formation energy ingredient
            Args:
                inputtxt <str>: Input text file, containing:
                    dfe_label1=perfect_label defected_label
                    dfe_label2=perfect_label defected_label
                    dfe_label3=perfect_label defected_label etc.
                    bandgap_lda_or_gga=<float>
                    bandgap_hse_or_expt=<float>
                plot_threshold <float>: Plotting threshold value
        """
        
        self.ingdir = os.getcwd() #Should be run in current directory
        self.recdir = os.path.dirname(self.ingdir)
        self.plot_threshold = plot_threshold
        ipparser = InputParser(inputfile=os.path.join(self.recdir,'input.inp'))
        self.input_options = ipparser.parse()
        self.dirs = dict()
        self.bandgap_lda_or_gga = 0.0
        self.bandgap_hse_or_expt = 0.0
        self.read_input_file(inputtxt)
        self.e_defects = dict()
        #pm = PickleManager(self.directory + '/input_options.pickle')
        #self.input_options = pm.load_variable()
        #from MAST.recipe.recipesetup import RecipeSetup
        #self.recipe_setup = RecipeSetup(inputOptions=self.input_options,recipeFile=os.path.join(self.recdir,'personal_recipe.txt'),workingDirectory=self.directory)
        #self.recipe_plan = self.recipe_setup.start()
        #PickleManager(self.directory + '/mast.pickle').load_variable()

        #self.final_ingredients = list()
        #ingredientlist = self.recipe_plan.ingredients.keys()
        #updatekeys = self.recipe_plan.update_methods.keys()
        #for ingredient in ingredientlist:
        #    if not (ingredient in updatekeys): # has no children
        #        self.final_ingredients.append(ingredient)
        #    else:
        #        if len(self.recipe_plan.update_methods[ingredient]) == 0:
        #            self.final_ingredients.append(ingredient)
        #print self.final_ingredients
        
    def read_input_file(self, inputtxt=""):
        """Read the input file.
            Args:
                inputtxt <str>: Input file name (shortname)
            Returns:
                Sets self.dirs with perfect and defected directories
                Sets self.bandgap_lda_or_gga
                Sets self.bandgap_hse_or_expt
        """
        ifpath = os.path.join(self.ingdir, inputtxt)
        ifile = MASTFile(ifpath)
        for dline in ifile.data:
            dline = dline.strip()
            if dline[0] == "#":
                continue
            if len(dline) == 0:
                continue
            dsplit = dline.split("=")
            dkey = dsplit[0]
            dval = dsplit[1]
            if dkey[0:4] == "dfe_":
                mylabel = dkey[4:]
                myperfect = dval.split()[0]
                mydefected = dval.split()[1]
                self.dirs[mylabel] = dict()
                self.dirs[mylabel]['perfect']=myperfect
                self.dirs[mylabel]['defected']=mydefected
            elif dkey == "bandgap_lda_or_gga"
                self.bandgap_lda_or_gga = float(dval)
            elif dkey == "bandgap_hse_or_expt"
                self.bandgap_hse_or_expt = float(dval)
            else:
                pass
        return

    def _calculate_defect_formation_energies(self):
        try:
            perf_dir = [ingredient for ingredient in self.final_ingredients if ('perfect' in ingredient)][0]
        except IndexError:
            raise MASTError(self.__class__.__name__, "A perfect final directory (has no children and has the word 'perfect' in its name) could not be found. Check recipe %s for a perfect directory." % self.directory)

        def_dir = [ingredient for ingredient in self.final_ingredients if 'perfect' not in ingredient] 
        #for whatever in sorted(def_dir):
            #print whatever

        defects = self.input_options.get_item('defects', 'defects')
        chempot = self.input_options.get_item('chemical_potentials')

        perf_meta = Metadata(metafile='%s/%s/metadata.txt' % (self.directory, perf_dir))
        e_perf = float(perf_meta.read_data('energy'))
        efermi = self.get_fermi_energy(perf_dir)
        struct_perf = self.get_structure(perf_dir)

        perf_species = dict()
        for site in struct_perf:
            if (str(site.specie) not in perf_species):
                perf_species[str(site.specie)] = 1
            else:
                perf_species[str(site.specie)] += 1

        #print 'Perfect composition: ', perf_species

        # First loop through the conditions for the defects
        for conditions, potentials in chempot.items():
            print 'Conditions:  %s' % conditions.upper()
            self.e_defects[conditions] = dict()

            # Loop through each defect
            for ddir in sorted(def_dir):
                def_meta = Metadata(metafile='%s/%s/metadata.txt' % (self.directory, ddir))
                #print def_meta
                label = def_meta.read_data('defect_label')
                charge = int(def_meta.read_data('charge'))
                energy = float(def_meta.read_data('energy'))
                structure = self.get_structure(ddir)

                if (label not in self.e_defects[conditions]):
                    self.e_defects[conditions][label] = list()

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
                #print '%-15s%-5i%12.5f%12.5f%12.5f%12.5f' % (label.split('_')[1], charge, energy, e_perf, efermi, alignment)
                print 'DFE = %f' % e_def
                self.e_defects[conditions][label].append( (charge, e_def) )

        #print e_defects
        #return e_defects

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
        elif ('OUTCAR' in os.listdir(abspath)):
            return Outcar('%s/OUTCAR' % abspath).efermi

    def get_structure(self, directory):
        """Returns the final structure from an optimization"""

        abspath = '%s/%s/' % (self.directory, directory)

        if ('vasprun.xml' in os.listdir(abspath)):
            return Vasprun('%s/vasprun.xml' % abspath).final_structure
        elif ('CONTCAR' in os.listdir(abspath)):
            return pmg.read_structure('%s/CONTCAR' % abspath)

    def get_potential_alignment(self, perf_dir, def_dir):
        """Returns the potential alignment correction used in charge defects"""

        abs_path_perf = '%s/%s/' % (self.directory, perf_dir)
        abs_path_def = '%s/%s/' % (self.directory, def_dir)

        pa = PotentialAlignment()

        if ('OUTCAR' in os.listdir(abs_path_perf)):
            perfect_info = pa.read_outcar('%s/%s' % (abs_path_perf, 'OUTCAR'))
            defect_info = pa.read_outcar('%s/%s' % (abs_path_def, 'OUTCAR'))

            return pa.get_potential_alignment(perfect_info, defect_info)

    def get_defect_formation_energies(self):
        if not self.e_defects:
            self._calculate_defect_formation_energies()

        return self.e_defects

    @property
    def defect_formation_energies(self):
        return self.get_defect_formation_energies()

    @property
    def dfe(self):
        return self.get_defect_formation_energies()

    def print_table(self):
        if not self.e_defects:
            self._calculate_defect_formation_energies()

        myfile = MASTFile()
        for conditions, defects in self.e_defects.items():
            myfile.data.append('\n\nDefect formation energies for %s conditions.\n' % conditions.upper())
            myfile.data.append('%-20s | %10s\n' % ('Defect', 'DFE'))
            myfile.data.append('---------------------------------\n')
            for defect, energies in defects.items():
                for energy in energies:
                    myfile.data.append('%-14s(q=%2i) | %8.4f\n' % (defect, energy[0], energy[1]))
                myfile.data.append(str()) # Add a blank line here
            myfile.data.append('---------------------------------\n')
        myfile.to_file(os.path.join(os.getcwd(),"dfe.txt"))
        for line in myfile.data:
            print line.strip()

