##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import sys, os

import pymatgen as pmg
from pymatgen.io.vasp import Vasprun
from MAST.utility import MASTError
from MAST.utility import MASTFile
from MAST.utility.defect_formation_energy.defectformationenergy import DefectFormationEnergy
from MAST.utility import loggerutils
from MAST.utility import fileutil
import logging
from MAST.ingredients.checker import VaspChecker
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
            self.e_defects <dict of list of tuples>: 
                Defect formation energy dictionary, by charge
                self.e_defects[conditions][label] =
                    [(charge <float>, dfe <float>),
                     (charge <float>, dfe <float), etc.]
            self.logger <logging logger>
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
        self.logger = loggerutils.get_mast_logger("DFE ingredient")
        self.ingdir = os.getcwd() #Should be run in current directory
        self.recdir = os.path.dirname(self.ingdir)
        self.plot_threshold = plot_threshold
        ipdir = os.path.join(self.recdir,'input.inp')
        os.chdir(self.recdir) #change to recipe directory
        ipparser = InputParser(inputfile=ipdir)
        self.input_options = ipparser.parse()
        os.chdir(self.ingdir) #change back to ingredient directory
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
            dsplit = dline.split("=",1)
            dkey = dsplit[0]
            dval = dsplit[1]
            if dkey[0:4] == "dfe_":
                mylabel = dkey[4:]
                myperfect = dval.split()[0]
                mydefected = dval.split()[1]
                self.dirs[mylabel] = dict()
                self.dirs[mylabel]['perfect']=myperfect
                self.dirs[mylabel]['defected']=mydefected
            elif dkey == "bandgap_lda_or_gga":
                self.bandgap_lda_or_gga = float(dval)
            elif dkey == "bandgap_hse_or_expt":
                self.bandgap_hse_or_expt = float(dval)
            else:
                pass
        return
    def calculate_defect_formation_energies(self):
        """Calculate all defect formation energies"""
        chempot = self.input_options.get_item('chemical_potentials')
        for conditions, potentials in chempot.items():
            print 'Conditions:  %s' % conditions.upper()
            self.e_defects[conditions] = dict()
            for mylabel in self.dirs.keys():
                self.calculate_single_defect_formation_energy(mylabel, conditions, potentials)


    def calculate_single_defect_formation_energy(self, mylabel, conditions, potentials):
        """Calculate one defect formation energy.
            Args:
                mylabel <str>: Defect label
                conditions <str>: Environment condition, e.g. "As-rich"
                potentials <dict of float>: Dictionary of chemical potentials,
                    e.g.:
                    potentials['As']=-6.0383
                    potentials['Ga']=-3.6080
                    potentials['Bi']=-4.5650
            Returns:
                <float>: Defect formation energy
        """
        pdir = self.dirs[mylabel]['perfect']
        ddir = self.dirs[mylabel]['defected']
        #Get perfect data
        perf_meta = Metadata(metafile='%s/%s/metadata.txt' % (self.recdir, pdir))
        e_perf = self.get_total_energy(pdir)
        #e_perf = float(perf_meta.read_data('energy').split(',')[-1]) #Sometimes the energy is listed more than once.
        efermi = self.get_fermi_energy(pdir)
        struct_perf = self.get_structure(pdir)
        perf_species = dict()
        for site in struct_perf:
            if (str(site.specie) not in perf_species):
                perf_species[str(site.specie)] = 1
            else:
                perf_species[str(site.specie)] += 1

        #print 'Perfect composition: ', perf_species

        #Get defected data
        def_meta = Metadata(metafile='%s/%s/metadata.txt' % (self.recdir, ddir))
        #print def_meta
        label = def_meta.read_data('defect_label')
        charge = int(def_meta.read_data('charge'))
        energy = self.get_total_energy(ddir)
        #energy = float(def_meta.read_data('energy').split(',')[-1])
        structure = self.get_structure(ddir)
        if (label not in self.e_defects[conditions]):
            self.e_defects[conditions][label] = list()
        print 'Calculating DFEs for defect %s with charge %3i.' % (label, charge)
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
        alignment = self.get_potential_alignment(pdir, ddir)

        # Calculate the base DFE energy
        e_def = energy - e_perf # E_defect - E_perf
        print "TTM DEBUG: e_def: ", e_def
        for specie, number in struct_diff.items():
            mu = potentials[str(specie)]
            #print str(specie), mu, number
            e_def -= (number * mu)
            print "TTM DEBUG: number: ", number
            print "TTM DEBUG: mu: ", mu
            print "TTM DEBUG: e_def -= number*mu: ", e_def
        print "TTM DEBUG: charge: ", charge
        print "TTM DEBUG: alignment: ", alignment
        print "TTM DEBUG: efermi: ", efermi
        e_def += charge * (efermi + alignment) # Add in the shift here!
        print "TTM DEBUG: e_def += charge*(efermi + alignment): ", e_def
        #print '%-15s%-5i%12.5f%12.5f%12.5f%12.5f' % (label.split('_')[1], charge, energy, e_perf, efermi, alignment)
        print 'DFE = %f' % e_def
        self.e_defects[conditions][label].append( (charge, e_def) )
        return

    def get_total_energy(self, directory):
        """Returns the total energy from a directory"""

        abspath = '%s/%s/' % (self.recdir, directory)
        mychecker = VaspChecker(name=abspath)
        return mychecker.get_energy_from_energy_file()

    def get_fermi_energy(self, directory):
        """Returns the Fermi energy from a directory"""
        abspath = '%s/%s/' % (self.recdir, directory)

        if ('OUTCAR' in os.listdir(abspath)):
            grepfermi = fileutil.grepme("%s/OUTCAR" % abspath, "E-fermi")
            lastfermi=grepfermi[-1]
            fermi = float(lastfermi.split()[2])
            print "TTM DEBUG LAST FERMI: ", fermi
            return fermi
        elif ('vasprun.xml' in os.listdir(abspath)):
            return Vasprun('%s/vasprun.xml' % abspath).efermi

    def get_structure(self, directory):
        """Returns the final structure from an optimization"""

        abspath = '%s/%s/' % (self.recdir, directory)

        if ('CONTCAR' in os.listdir(abspath)):
            return pmg.core.structure.IStructure.from_file('%s/CONTCAR' % abspath)
        elif ('vasprun.xml' in os.listdir(abspath)):
            return Vasprun('%s/vasprun.xml' % abspath).final_structure

    def get_potential_alignment(self, perf_dir, def_dir):
        """Returns the potential alignment correction used in charge defects"""

        abs_path_perf = '%s/%s/' % (self.recdir, perf_dir)
        abs_path_def = '%s/%s/' % (self.recdir, def_dir)

        pa = PotentialAlignment()

        if ('OUTCAR' in os.listdir(abs_path_perf)):
            perfect_info = pa.read_outcar('%s/%s' % (abs_path_perf, 'OUTCAR'))
            defect_info = pa.read_outcar('%s/%s' % (abs_path_def, 'OUTCAR'))
            print "TTM DEBUG: Perfect PA info: ", perfect_info
            print "TTM DEBUG: Defected PA info: ", defect_info
            pa_output = pa.get_potential_alignment(perfect_info, defect_info)
            print "TTM DEBUG: PA output: ", defect_info
            return pa_output

    def get_defect_formation_energies(self):
        if not self.e_defects:
            self.calculate_defect_formation_energies()

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
            self.logger.info(line.strip())


