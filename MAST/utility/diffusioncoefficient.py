#!/usr/bin/env python
import os
import sys
import pymatgen
import numpy as np
from MAST.utility import MASTError
from MAST.utility import MASTFile
from MAST.utility import PickleManager
from MAST.utility import dirutil
kboltz = 8.6173325E-5
stock_v=1e13
stock_S_v=2.07e-4 #Shewmon, Sv/R approx 2.4?? =2.07e-4
stock_S_m=2.07e-4

class DiffusionCoefficient():
    """Dilute solute vacancy diffusion coefficient calculations.
        DOES NOT SUPPORT CHARGED SUPERCELLS until DefectFormationEnergy does
        Args:
            self.directory <str>: directory of the recipe to calculate
            self.freqmodel <int>: frequency model   
                                5 : five-frequency model
                                1 : one frequency; pure system
            self.freqdict <dict of str>: frequency-to-label dictionary:
                        freqdict['w0']=['vac1-vac2']
                        freqdict['w4']=['vac10-vac9']
                            If a label is backwards from what is in the 
                            neb ingredient's actual name, then
                            the backwards hop will be taken.
                            Initial state for phonons will be the first
                            label in the label pair as listed in this dict.
            self.host <str>: string with the host element's symbol (e.g. "Fe")
            self.solute <str>: string with the solute element's symbol
            self.input_options <InputOptions>: input options of the recipe
            self.dirlist <list of str>: list of directories underneath 
                                        self.directory
            self.hopdict <dict of float>: freq-labeled dict of NEB barriers
                                        (right now this is ENERGIES, NOT
                                        ENTHALPIES)
            self.attemptdict <dict of float>: freq-labeled dict of attempt freqs
            self.entromdict <dict of float>: freq-labeled dict of migration entropies
            self.formdict <dict of float>: freq-labeled dict of formation
            self.tempdict <dict of float>: temperature list of diffusion coeffs
    """
    def __init__(self, directory, tempstart=73, tempend=1273, tempstep=100, freqmodel=0, freqdict=None, verbose=0):
        """
            Args:
                directory <str>: directory to start in
                tempstart <float>: Temp in K to start range (default 73)
                tempend <float>: Temp in K to end range (default 1273)
                tempstep <float>: step in K for temp. range (default 100)
                freqmodel <int>: Integer for frequency model
                freqdict <dict>: Dictionary for frequency-to-hop-labels
                verbose <int>: 0 for not verbose (default), 1 for verbose
        """
        self.directory = directory
        dirlist = os.listdir(self.directory)
        dirlist.sort()
        self.dirlist = dirlist
        self.freqmodel=None
        self.freqdict=None
        self.host=None
        self.solute=None
        self.input_options=None
        self.hopdict=None
        self.attemptdict=None
        self.entrodict=None
        self.formdict=None
        self.tempdict=None
        self.verbose=verbose
        if self.verbose == 1:
            print "verbose mode is on"
        try:
            pm = PickleManager(os.path.join(self.directory,
                                        'input_options.pickle'))
            self.input_options = pm.load_variable()
        except IOError:
            pass
        self.set_up_temp_dict(tempstart, tempend, tempstep)
        self.set_freqmodel(freqmodel)
        self.set_frequency_dict(freqdict)
        self.set_host_and_solute()

        return

    def set_up_temp_dict(self, tempstart=73, tempend=1273, tempstep=100):
        """Set up temperature dict.
            Args:
                tempstart <float>: starting temp in K (default 73 K)
                tempend <float>: ending temp in K (default 1273 K)
                tempstep <float>: step in K (default 100 K)
        """
        if tempend < tempstart:
            tempend = tempstart #Do not allow bad temp range.
        if tempstart <= 0:
            tempstart = 0.000001 #Do not allow negative or zero temperature.
        tidx = tempstart
        tdict=dict()
        tct=0
        maxidx=10000 #arbitrary stop at 10000 in case of errors
        while (tidx <= tempend) and (tct < maxidx):
            tdict[tidx] = 0.
            tidx = tidx + tempstep
            tct = tct + 1
        self.tempdict = tdict
        if self.verbose == 1:
            print "Temp dict: ", self.tempdict
        return self.tempdict

    def set_freqmodel(self, freqmodel=0):        
        """Set the frequency model to use.
            Args:
                freqmodel <int>:
                                1 : pure vacancy
                                5 : five-frequency model
                                0 : set from input options if possible.
        """
        trymodel=None
        if (freqmodel == 0) and not (self.input_options == None):
            if "freqmodel" in self.input_options['program_keys'].keys():
                trymodel = self.input_options['program_keys']['freqmodel']
        if trymodel == None: #was not set in input options
            try:
                trymodel=int(freqmodel)
            except (TypeError, ValueError):
                raise MASTError(self.__class__.__name__, "Bad type for frequency model integer.")
        if not (trymodel in [1, 5]):
            raise MASTError(self.__class__.__name__, "%s is not a supported frequency model." % str(trymodel))
        self.freqmodel = trymodel
        if self.verbose == 1:
            print "Freq model: ", self.freqmodel
        return self.freqmodel
            
    def set_frequency_dict(self, freqdict=None):
        """Set up the frequency dictionary.
            Args:
                freqdict <dict of str>: Frequency dictionary. 
                    This is a frequency-to-label dictionary with the format:
                        freqdict['w0']=['vac1-vac2']
                        freqdict['w4']=['vac10-vac9']
                            If a label is backwards from what is in the 
                            neb ingredient's actual name, then
                            the backwards hop will be taken.
                            Initial state for phonons will be the first
                            label in the label pair as listed in this dict.
                    If None, try to set from input options
        """
        trydict=None
        if (freqdict == None) and not (self.input_options == None):
            if "freqdict" in self.input_options['program_keys'].keys():
                trydict = dict(self.input_options['program_keys']['freqdict'])
        if trydict==None: #was not set in input options
            trydict = dict(freqdict)
        if trydict==None:
            raise MASTError(self.__class__.__name__, "No frequency dictionary.")
        for (key, value) in trydict.iteritems(): #validate
            if not 'w' in key:
                raise MASTError(self.__class__.__name__, "Frequency dict does not have the correct format, e.g. freqdict['w0']=['vac1-vac2']")
            if not '-' in value:
                raise MASTError(self.__class__.__name__, "Frequency dict does not have the correct format, e.g. freqdict['w0']=['vac1-vac2']")
        self.freqdict = trydict
        if self.verbose == 1:
            print "Freq dict: ", self.freqdict
        return self.freqdict


    def set_host_and_solute(self):
        """Find the host and the solute based on relative numbers in POSCAR,
            for a BINARY, DILUTE system. (One host element type, 
            one solute element type.)
        """
        label=""
        if self.freqmodel == 1:
            if not 'w0' in self.freqdict.keys():
                raise MASTError(self.__class__.__name__, "Need a label for frequency w0 for one-frequency model.")
            label = self.freqdict['w0']
        else:
            if not 'w2' in self.freqdict.keys():
                raise MASTError(self.__class__.__name__, "Not enough frequencies to determine solute and host.")
            label = self.freqdict['w2']       
        [startlabel, tstlabel]=self.find_start_and_tst_labels(label)
        mydir = self.get_labeled_dir(startlabel, "stat","","neb")
        pospath = os.path.join(mydir, "POSCAR")
        if not os.path.isfile(pospath):
            raise MASTError(self.__class__.__name__, "No POSCAR in %s." % pospath)
        mypos = pymatgen.io.vaspio.Poscar.from_file(pospath)
        sitesym = mypos.site_symbols
        natoms = mypos.natoms
        hostidx=natoms.index(max(natoms))
        self.host=sitesym[hostidx]
        if len(natoms) > 1:
            solidx=natoms.index(min(natoms))
            self.solute=sitesym[solidx]
        if self.verbose == 1:
            print "host: ", self.host
            print "solute: ", self.solute
        return 

    def find_start_and_tst_labels(self, label, app1="neb", app2="stat"):
        """Finds the start label and transition state labels for ingredients.
            Args:
                label <str>: neb label, like "vac1-vac2" 
                app1 <str>: additional string which should be in the label
                app2 <str>: additional string which should be in the label
            Returns:
                [startlabel, tstlabel]
                startlabel <str>: label for start directory (e.g. "vac1" or "vac2")
                tstlabel <str>: actual directory label for neb (e.g. "vac1-vac2" or "vac2-vac1")
        """
        startlabel=""
        tstlabel=""
        lsplit = label.split('-')
        backlab = '-'.join([lsplit[1],lsplit[0]])
        for myname in self.dirlist:
            if (backlab in myname) and (app1 in myname) and (app2 in myname):
                tstlabel = backlab
                startlabel = label.split('-')[0]
                break
            elif (label in myname) and (app1 in myname) and (app2 in myname):
                tstlabel = label
                startlabel = label.split('-')[0]
                break
        if tstlabel == "":
            raise MASTError(self.__class__.__name__, 
                "No directory found for hop label %s." % label)
        if self.verbose == 1:
            print "Transition start and transition state: ", [startlabel, tstlabel]
        return [startlabel, tstlabel]

    def get_labeled_dir(self, label, append="stat", append2="", ignore=""):
        """Get a labeled directory.
            Args:
                label <str>: label for the directory, like "vac1"
                append <str>: extra piece which should appear in the 
                                directory name
                append2 <str>: additional extra piece which should appear
                ignore <str>: piece to ignore
            Returns:
                mypath <str>: The FIRST matching directory found, from a 
                                sorted list. Returns a full path.
        """
        dirlist = os.listdir(self.directory)
        dirlist.sort()
        mydir=""
        mypath=""
        for mydir in dirlist:
            if (label in mydir) and (append in mydir) and (append2 in mydir):
                if not (ignore == "") and (ignore in mydir):
                    pass
                else:
                    mypath=os.path.join(self.directory, mydir)
                    break
        if self.verbose == 1:
            print "mypath: ", mypath
        if mypath=="":
            raise MASTError(self.__class__.__name__, "No path found for label %s with append %s %s" % (label,append,append2))
        return mypath

    def neb_barrier(self, label):
        """Get the NEB barrier for one label.
            If the label is reversed from the folder name, the "start"
            directory will be the endpoint whose label comes first.
            The transition state is the maximum energy.
            Performs NO CHECKS for the state of the energy contour.
            Args:
                label <str>: label for the NEB folder, like "vac1-vac2"
            Returns:
                nebenergy <float>: transition state energy-starting state energy
        """
        [startlabel, tstlabel] = self.find_start_and_tst_labels(label)
        startenergy = self.get_labeled_energy(startlabel,"stat","","neb")
        maxenergy = self.get_max_energy(tstlabel)
        if self.verbose == 1:
            print "Start energy: ", startenergy
            print "Max energy: ", maxenergy
        if (startenergy > maxenergy):
            print "Attention! Maximum energy is from an endpoint, not an image!"
        nebenergy = maxenergy - startenergy
        if self.verbose == 1:
            print "NEB energy: ", nebenergy
        return nebenergy

    def get_labeled_energy(self, label, append="stat", append2="", ignore=""):
        """Get an energy from a directory.
            Args:
                label <str>: label for directory
                append <str>: additional tag for the directory
                append2 <str>: additional tag for the directory
                ignore <str>: string to ignore
        """
        mypath = self.get_labeled_dir(label, append, append2, ignore)
        myenergy=0
        myenergy = self.get_total_energy(mypath)
        if self.verbose == 1:
            print "Myenergy: ", myenergy
        return myenergy

    def get_total_energy(self, directory): #synchronize this call with Glen's 
        """Returns the total energy from a directory.
            Args:
                directory <str>: directory (full path)
            Returns:
                <float>: Total energy (E0 from VASP)
            SHOULD BE REPLACED WITH PROGRAM-INDEPENDENT SWITCH.
        """
        if 'OSZICAR' in os.listdir(directory):
            ozfile = MASTFile(os.path.join(directory, "OSZICAR"))
            ezeroenergy = float(ozfile.get_segment_from_last_line_match("E0",
                    "E0= ","  d E"))
            if self.verbose == 1:
                print "E0 energy: ", ezeroenergy
            return ezeroenergy

    def get_max_energy(self, tstlabel, app1="neb", app2="stat"):
        """Get maximum image energy from an NEB directory.
            Args:
                tstlabel <str>: label for the NEB directory
                app1 <str>: label which must be in the name.
                app2 <str>: label which  must be in the name
            Returns:
                <float>: Maximum energy among the images.
        """
        mypath = self.get_labeled_dir(tstlabel, app1, app2)
        imdirs = os.listdir(mypath)
        imdirs.sort()
        mydir=""
        imenergs=dict()
        if self.verbose == 1:
            print "imdirs: ", imdirs
        for mydir in imdirs:
            fulldir = os.path.join(mypath, mydir)
            if os.path.isdir(fulldir):
                if self.verbose == 1:
                    print "using dir: ", fulldir
                imenergs[self.get_total_energy(fulldir)] = mydir
        maxenerg = max(imenergs)
        #maxdir = imenergs[maxenerg]
        if self.verbose == 1:
            print "max energy: ", maxenerg
        return maxenerg
   

    def get_hopdict(self):
        """Get dictionary of hops."""
        freqs = self.freqdict.keys()
        hopdict=dict()
        for freq in freqs:
            neblabel = self.freqdict[freq]
            hopdict[freq] = self.neb_barrier(neblabel)
        self.hopdict = hopdict
        if self.verbose == 1:
            print "Hopdict: ", self.hopdict
        return self.hopdict

    def get_formdict(self):
        """Get formation energy dictionary.
        """
        formdict=dict()
        #pure_def = #TTM DEBUG NEED TO SWITCH TO USE GLEN'S FUNCTION
        for freq in self.freqdict.keys():
            formdict[freq] = self.get_defect_energy_simple(self.freqdict[freq])
        self.formdict = formdict
        if self.verbose == 1:
            print "Formdict: ", self.formdict
        return self.formdict

    def get_defect_energy_simple(self, label, append="stat", perf="perf"):
        """Need to replace with Glen's function.
            Get simplified defect formation energy. NO CHARGED CELLS.
            Args:
                label <str>: label for hop (starting label will be determined
                            and then compared against the pure).
                append <str>: append for name 
                perf <str>: perfect-cell designation
            Return:
                <float>: defect formation energy
        """
        [startlabel, tstlabel]=self.find_start_and_tst_labels(label)
        defpath=self.get_labeled_dir(startlabel, append, "", "neb")
        perfpath=self.get_labeled_dir(perf, append)
        defectedenergy = self.get_total_energy(defpath)
        perfectenergy = self.get_total_energy(perfpath)
        if defectedenergy==None or perfectenergy==None:
            raise MASTError(self.__class__.__name__,"No defected energy found.")
        delta = defectedenergy - perfectenergy
        refpure=0 # Need to get this from a calc??
        if refpure==0:
            mypos=pymatgen.io.vaspio.Poscar.from_file(os.path.join(perfpath, "POSCAR"))
            totatoms=sum(mypos.natoms)
            refpure=perfectenergy/totatoms
            #myfile=MASTFile()
            #myfile.data.append("WARNING: Cohesive energy has not been properly calculated in %s because a pseudopotential reference has not been taken into account." % perfpath)
            #myfile.to_file(os.path.join(self.directory,"WARNING"))
        eform = delta + refpure #Add on reference atom
        if self.verbose == 1:
            print "Eform: ", eform
        return eform

    def get_attemptdict(self):
        """Get approximate attempt vibrational frequency dictionary, 
            using Adams approximations:
            v_0=v_1=v_3=v_4 and (v2/v0)=sqrt(mass0*Tmeltpure2/(mass2*Tmeltpure0)
            Temperatures are in KELVIN.
            Masses are in AMU.
        """
        vibdict=dict()
        freqs = self.freqdict.keys()
        freqs.sort()
        for freq in freqs:
            if not (freq == 'w2'):
                vibdict[freq] = stock_v
            else:
                hostelem = pymatgen.core.periodic_table.Element(self.host)
                solelem = pymatgen.core.periodic_table.Element(self.solute)
                masshost = hostelem.atomic_mass
                masssolute = solelem.atomic_mass
                tmeltpurehost = self.parse_melting_point(hostelem)
                tmeltpuresolute = self.parse_melting_point(solelem)
                vibdict[freq] = vibdict['w0'] * np.sqrt(masshost*tmeltpuresolute/(masssolute*tmeltpurehost))
        self.attemptdict=vibdict
        if self.verbose == 1:
            print "attemptdict: ", self.attemptdict
        return self.attemptdict
    

    def get_entromigdict_approx(self):
        """Get approximate entropy.
            Use Adams approximations for now, Wert-Zener, 
            where S_0=S_1=S_3=S_4 and
            S_2 - S_0 = beta(E_2 - E_1)/Tmelt, pure,
            where S are activation entropies and
            E are activation energies and beta is assumed to be 0.4

            Dobson et al, S/k approx 10?            
            Shewmon, Sm approx beta*Q/Tm, Sv/R = 2.4
            http://www.tf.uni-kiel.de/matwis/amat/
                def_en/kap_3/backbone/r3_1_1.html, Sform ~ Smig ~ 1k?         
            LATER, REPLACE WITH EXTRACTION of "THERMO" from PHON by Dario Alfe
        """
        entromdict=dict()
        beta=0.4
        hostelem = pymatgen.core.periodic_table.Element(self.host)
        tmeltpurehost = self.parse_melting_point(hostelem)
        freqs = self.freqdict.keys()
        freqs.sort() #do w0 before w2
        for freq in freqs:
            if not (freq == 'w2'):
                entromdict[freq] = stock_S_m
            else:
                eact2 = self.hopdict['w2']
                eact1 = self.hopdict['w1']
                delta_act_s = beta*(eact2-eact1)/tmeltpurehost
                entromdict['w2'] = entromdict['w0'] + delta_act_s
        self.entromdict = entromdict
        if self.verbose == 1:
            print "Entro mig dict: ", self.entromdict
        return self.entromdict

    def get_entrodict_from_phonons(self, temp):
        """Get entropy dictionary from phonon calculations, integrated by
            PHON.
            Args:
                temp <float>: Temperature in K
        """
        entromdict=dict() #migration entropy
        entrovdict=dict() #vacancy formation entropy
        entrodict=dict()
        freqs = self.freqdict.keys()
        freqs.sort() #do w0 before w2
        bulkfind=dirutil.search_for_metadata_file("phononlabel=perfect,ingredient type=PhonParse", self.directory)
        if len(bulkfind) > 1:
            raise MASTError(self.__class__.__name__, "More than one bulk phonon parse found.")
        bulkfile = MASTFile(os.path.join(bulkfind[0] + "THERMO"))
        entropybulk=0
        for myline in bulkfile.data:
            mlist = myline.strip().split()
            if float(mlist[0]) == temp:
                entropybulk = float(mlist[4]) #in kB/cell
                break
        for freq in freqs:
            [startlabel,tstlabel]=self.find_start_and_tst_labels(self.freqdict[freq])
            initfind=dirutil.search_for_metadata_file("phononlabel=" + startlabel + ",ingredient type=PhonParse", self.directory)
            tstfind=dirutil.search_for_metadata_file("phononlabel=" + tstlabel + ",ingredient type=PhonParse", self.directory)
            if len(initfind) > 1:
                raise MASTError(self.__class__.__name__, "More than one initial phonon parse found for %s" % startlabel)
            if len(tstfind) > 1:
                raise MASTError(self.__class__.__name__, "More than one transition phonon parse found for %s" % tstlabel)
            initfile = MASTFile(os.path.join(os.path.dirname(initfind[0]),"THERMO"))
            tstfile = MASTFile(os.path.join(os.path.dirname(tstfind[0]),"THERMO"))
            entropyinit=0
            entropytst=0
            for myline in initfile.data:
                mlist = myline.strip().split()
                if float(mlist[0]) == temp:
                    entropyinit = float(mlist[4]) #in kB/cell
                    break
            for myline in tstfile.data:
                mlist = myline.strip().split()
                if float(mlist[0]) == temp:
                    entropytst = float(mlist[4]) #in kB/cell
                    break
            entromdict[freq]=(entropytst-entropyinit)
            entrovdict[freq]=(entropyinit-entropybulk)
            entrodict[freq]=entrovdict[freq]+entromdict[freq]
        self.entrodict = entrodict
        if self.verbose == 1:
            print "Entrodict: ", self.entrodict
        return self.entrodict

    def get_entromigdict_from_phonons(self, temp):
        """Get entropy migration dictionary from phonon 
            calculations, integrated by PHON.
            Args:
                temp <float>: Temperature in K
        """
        entromdict=dict() #migration entropy
        freqs = self.freqdict.keys()
        freqs.sort() #do w0 before w2
        for freq in freqs:
            [startlabel,tstlabel]=self.find_start_and_tst_labels(self.freqdict[freq])
            initfind=dirutil.search_for_metadata_file("phononlabel=" + startlabel + ",ingredient type=PhonParse", self.directory)
            tstfind=dirutil.search_for_metadata_file("phononlabel=" + tstlabel + ",ingredient type=PhonParse", self.directory)
            if self.verbose == 1:
                print "initfind: ", initfind
                print "tstfind: ", tstfind
            if len(initfind) > 1:
                raise MASTError(self.__class__.__name__, "More than one initial phonon parse found for %s" % startlabel)
            if len(tstfind) > 1:
                raise MASTError(self.__class__.__name__, "More than one transition phonon parse found for %s" % tstlabel)
            initfile = MASTFile(os.path.join(os.path.dirname(initfind[0]),"THERMO"))
            tstfile = MASTFile(os.path.join(os.path.dirname(tstfind[0]),"THERMO"))
            entropyinit=0
            entropytst=0
            mfloat=0.0
            for myline in initfile.data:
                mlist = myline.strip().split()
                try: 
                    mfloat=float(mlist[0])
                    if mfloat == temp:
                        entropyinit = float(mlist[4]) #in kB/cell
                        break
                except (ValueError):
                    pass
            for myline in tstfile.data:
                mlist = myline.strip().split()
                try: 
                    mfloat=float(mlist[0])
                    if mfloat == temp:
                        entropytst = float(mlist[4]) #in kB/cell
                        break
                except (ValueError):
                    pass
            entromdict[freq]=(entropytst-entropyinit)
        self.entromdict = entromdict
        if self.verbose == 1:
            print "Last entry for entropytst: ", entropytst
            print "Last entry for entropyinit: ", entropyinit
            print "Entromigdict: ", self.entromdict
        return self.entromdict
    def parse_melting_point(self, elem):
        """Parse the melting point.
            Args:
                elem <Element>: pymatgen Element class object
            Return:
                <float>: melting point in Kelvin
        """
        meltstr = str(elem.melting_point) #str(<unicode>)
        meltval = meltstr.split()[0]
        meltval = float(meltval)
        if self.verbose == 1:
            print "Melting point: ", meltval
        return meltval


    def five_freq(self, temp):
        """Five-frequency model for FCC SYSTEMS ONLY.
            Args:
                temp <float>: Temperature in K
            Returns:
                [Dself <float>, Dsolute <float>]
                    Self-diffusion coefficient
                    Dilute solute vacancy diffusion coefficient (tracer?)
        Adams, J. B., Foiles, S. M. & Wolfer, W. G. 
        Self-diffusion and impurity diffusion of fcc metals 
        using the five-frequency model and the Embedded Atom Method. 
        Journal of Materials Research 4, 102.112 (2011)
        """
        self.get_hopdict()
        self.get_formdict()
        #self.get_entromigdict_approx()
        self.get_entromigdict_from_phonons(temp)
        self.get_attemptdict() #Adams approximations
        freqs = self.freqdict.keys()
        jumpfreqdict=dict()
        for freq in freqs:
            jumpfreqdict[freq] = self.estimate_jump_freq(self.attemptdict[freq],
                self.entromdict[freq], self.hopdict[freq], temp)
        jfw0 = jumpfreqdict['w0']
        jfw1 = jumpfreqdict['w1']
        jfw2 = jumpfreqdict['w2']
        jfw3 = jumpfreqdict['w3']
        jfw4 = jumpfreqdict['w4']
        f_0 = 0.7815
        myx = jfw4/jfw0
        bigf_num = 10*np.power(myx,4) + 180.5*np.power(myx,3) + 927*np.power(myx,2) + 1341*myx
        bigf_denom = 2*np.power(myx,4) + 40.2*np.power(myx,3) + 254*np.power(myx,2) + 597*myx + 435
        bigf_x = 1-(1.0/7.0)*(bigf_num/bigf_denom)
        f_2_num = 1+3.5*bigf_x *(jfw3/jfw1)
        f_2_denom = 1+(jfw2/jfw1) + 3.5*bigf_x*(jfw3/jfw1)
        f_2 = f_2_num / f_2_denom
        vacconc = self.get_vac_conc(temp)
        lattparam = self.get_latt_param()
        Dself = f_0*vacconc*np.power(lattparam,2)*jfw0
        Dsol_num = Dself * f_2 * jfw2 * jfw4 * jfw1
        Dsol_denom = f_0 * jfw1 * jfw0 * jfw3
        Dsolute = Dsol_num / Dsol_denom
        if self.verbose == 1:
            print "Dself, Dsolute: ", [Dself, Dsolute]
        return [Dself, Dsolute]

    def get_vac_conc_cell(self, perf="perf"):
        """Get vacancy concentration of SUPERCELL. SHOULD WE USE EQUILIBRIUM VACANCY CONC.
            INSTEAD OF CELL VACANCY CONC?
            Args:
                perf <str>: designation string for perfect cell
                            (default) "perf"
            Returns:
                <float> vacancy concentration (UNITLESS)
        """
        label = self.freqdict['w0']       
        [startlabel, tstlabel]=self.find_start_and_tst_labels(label)
        defpath=self.get_labeled_dir(startlabel)
        defpospath = os.path.join(defpath, "POSCAR")
        if not os.path.isfile(defpospath):
            raise MASTError(self.__class__.__name__, "No POSCAR in %s." % defpospath)
        defpos = pymatgen.io.vaspio.Poscar.from_file(defpospath)
        defnatoms = defpos.natoms
        perfpath=self.get_labeled_dir(perf)
        perfpospath = os.path.join(perfpath, "POSCAR")
        if not os.path.isfile(perfpospath):
            raise MASTError(self.__class__.__name__, "No POSCAR in %s." % perfpospath)
        perfpos = pymatgen.io.vaspio.Poscar.from_file(perfpospath)
        perfnatoms = perfpos.natoms
        vacconc = (perfnatoms - defnatoms) / perfnatoms #e.g. (32-31)/32
        if self.verbose == 1:
            print "Cell vac conc: ", vacconc
        return vacconc

    def get_vac_conc(self, temp):
        """Get vacancy concentration in equilibrium.
            Shewmon 2nd Edition p. 70# N_v_eq = exp(Sv/R)exp(-Hv/RT)
            Here we use N_v_eq = exp(Sv/k)exp(-Hv/kT) with Sv in eV/K, Hv in eV
            Args:
                temp <float>: temp in K
            Returns:
                <float> vacancy concentration (UNITLESS)
        """
        label = self.freqdict['w0']       
        defenergy=self.get_defect_energy_simple(label)
        nveq = np.exp(stock_S_v/kboltz)*np.exp(-1*defenergy/(kboltz*temp))
        if self.verbose == 1:
            print "Equilib vac conc: ", nveq
        return nveq

    def get_latt_param(self, perf="perf"):
        """Get lattice parameter
            Args:
                perf <str>: designation string for perfect cell
                            (default) "perf"
            Returns:
                <float>: lattice parameter in ***centimeters*** (not Angstroms).
        """
        perfpath=self.get_labeled_dir(perf)
        perfpospath = os.path.join(perfpath, "POSCAR")
        if not os.path.isfile(perfpospath):
            raise MASTError(self.__class__.__name__, "No POSCAR in %s." % perfpospath)
        perfpos = pymatgen.io.vaspio.Poscar.from_file(perfpospath)
        numatoms = sum(perfpos.natoms)
        divisor=0
        try:
            closediv = np.power(numatoms/4, 0.33333333)
        except (ValueError, TypeError):
            raise MASTError(self.__class__.__name__, "Number of atoms could not be determined..")
        if np.abs(np.round(closediv) - closediv) > 0.001:
            raise MASTError(self.__class__.__name__, "This is apparently not a cubic lattice-vectored FCC cell and the divisor cannot be determined.")
        divisor = np.round(closediv) 
        #if numatoms == 108:
        #    divisor = 3 #This is a 3x3x3 unit cell
        #elif numatoms == 32:
        #    divisor = 2 #This is a 2x2x2 unit cell
        #elif numatoms == 4: 
        #    divisor = 1 #This is a single 'cubic-vectored' atom
        perfstr = perfpos.structure
        mylatt = float(perfstr.lattice.a)/float(divisor)
        mylatt = mylatt*np.power(10.0,-8)
        if self.verbose == 1:
            print "Latt param in centimeters: ", mylatt
        return mylatt
    
    def estimate_jump_freq(self, vibfreq, migentropy, migenergy, temp):
        """Estimate jump frequency (Adams, Foiles, Wolfer)
            w_i = v_i exp(S_i/k) exp(-E_i/kT)
            Use only Smig, Emig in here for five-frequency.
            Evf, Svf are taken in "vacancy concentration" (see five_freq)
        """
        jumpfreq = vibfreq*np.exp(migentropy/kboltz)*np.exp(-1*migenergy/(kboltz*temp))
        if self.verbose == 1:
            print "Act energy: ", migenergy
            print "Jump freq: ", jumpfreq
        return jumpfreq

    def diffusion_coefficient(self):
        """Diffusion coefficient switcher and temperature looper."""
        tkeys = self.tempdict.keys()
        tkeys.sort()
        for tkey in tkeys:
            if self.verbose == 1:
                print "Temperature (K): %1i\n" % tkey
            if self.freqmodel == 1:
                self.tempdict[tkey] = self.one_freq(tkey)
            elif self.freqmodel == 5:
                self.tempdict[tkey] = self.five_freq(tkey)
            else:
                raise MASTError(self.__class__.__name__, "%s is not a supported frequency model." % int(self.freqmodel))
        if self.verbose == 1:
            print "temp dict of diff coeffs: ", self.tempdict
        return self.tempdict


    def one_freq(self, temp):
        """One-frequency (pure) model for FCC SYSTEMS ONLY.
            Args:
                temp <float>: Temperature in K
            Returns:
                Dself <float>
                    Self-vacancy-diffusion coefficient
        """
        self.get_hopdict()
        self.get_formdict()
        #self.get_entromigdict_approx()
        self.get_entromigdict_from_phonons(temp)
        self.get_attemptdict() #Adams approximations
        freqs = self.freqdict.keys()
        jumpfreqdict=dict()
        for freq in freqs:
            jumpfreqdict[freq] = self.estimate_jump_freq(self.attemptdict[freq],
                self.entromdict[freq],
                self.hopdict[freq], temp)
        jfw0 = jumpfreqdict['w0']
        #f_0 = 0.715
        vacconc = self.get_vac_conc(temp)
        lattparam = self.get_latt_param()
        Dself = vacconc*np.power(lattparam,2)*jfw0
        if self.verbose == 1:
            print "Dself: ", Dself
        return Dself

    def print_temp_dict(self):
        """Print the self.tempdict dictionary of coefficients."""
        print "-----------------------------------------------"
        if self.freqmodel > 1:
            print "Host: ", self.host
            print "Solute: ", self.solute
            print "%8s %20s %20s" % ("TEMP (K)","Dself cm2/sec","Dsolute cm2/sec")   
        else:
            print "Pure: ", self.host
            print "%8s %20s" % ("TEMP (K)","D cm2/sec")   
        tkeys = self.tempdict.keys()
        tkeys.sort()
        for (tkey) in tkeys:
            tval = self.tempdict[tkey]
            if self.freqmodel > 1:
                print "%8s %17.2e %17.2e" % (tkey, tval[0], tval[1])
            else:
                print "%8s %17.2e" % (tkey, tval)
        print "-----------------------------------------------"

        return
def main():
    """Run diffusion coefficient calculator from the prompt.
        Args (sys.argv):
            directory <str>: path to recipe
            tstart <float>: start temp
            tend <float>: end temp
            tstep <float>: temp step
            freqmodel <int>: frequency model
            freqdict <str>: string to be converted into a frequency dict, of
                format "w0=vac1-vac2,w1=vac1-vac9,w2=vac1-vac3,w3=vac4-vac5..."
    """
    print "Arguments: ", sys.argv
    try:
        directory=os.path.expandvars(sys.argv[1])
    except IndexError:
        print sys.argv
    tstart=73
    tend=1273
    tstep=100
    freqmodel=1
    try:
        trytstart=sys.argv[2]
        tstart=float(trytstart)
    except IndexError:
        pass
    try:
        trytend=sys.argv[3]
        tend=float(trytend)
    except IndexError:
        pass
    try:
        trytstep=sys.argv[4]
        tstep=float(trytstep)
    except IndexError:
        pass
    try:
        tryfreqmodel=sys.argv[5]
        freqmodel=int(tryfreqmodel)
    except IndexError:
        pass
    freqdict=None
    try:
        fdstr = sys.argv[6]
        freqdict=dict()
        fdlist=fdstr.strip().split(",")
        for fditem in fdlist:
            mykey = fditem.split("=")[0]
            myval = fditem.split("=")[1]
            freqdict[mykey]=myval
    except IndexError:
        pass
    verbose=0
    try:
        verbose=int(sys.argv[7])
    except IndexError:
        pass

    print 'Looking at diffusion coefficient for %s' % directory
    DC = DiffusionCoefficient(directory, tstart, tend, tstep, freqmodel, freqdict, verbose)
    DC.diffusion_coefficient()
    DC.print_temp_dict()

if __name__ == '__main__':
    main()

#
def no_routine():
    """
            # Process NEB runs. The energy and magnetic moments for each image are listed.
            # Some useful values are calculated.

            # Number of NEB images
            image_list = []
            for i in xrange(0,9):
                append = "0"+str(i)
                newpath = os.path.join(fullpath,append)
                if os.path.exists(newpath):
                    image_list.append(newpath)
            d["num_images"] = len(image_list)

            # Image energies and magnetic moments for specific folders
            list_image_energies = []
            list_image_mags = []
            for i in xrange(0,len(image_list)):
                append = "0"+str(i)
                oszicar = os.path.join(fullpath,append,"OSZICAR")
                if not os.path.isfile(oszicar):
                    return None
                val_energy = util2.getEnergy(oszicar)
                val_mag = util2.getMag(oszicar)
                d["E_"+append]= val_energy
                d["mag_"+append]= val_mag
                list_image_energies.append(val_energy)
                list_image_mags.append(val_mag)

            # List of image energies and magnetic moments in order
            image_energies = ' '.join(map(str,list_image_energies))
            image_mags = ' '.join(map(str,list_image_mags))
            d["image_energies"] = image_energies
            d["image_mags"] = image_mags

            # An simple way to visualize relative image energies and magnetic moments
            energy_contour = "-x-"
            if len(image_list)==0:
                return None
            for i in xrange(1,len(image_list)):
                if(list_image_energies[i]>list_image_energies[i-1]):
                    energy_contour += "/-x-"
                elif list_image_energies[i]<list_image_energies[i-1]:
                    energy_contour += "\\-x-"
                else:
                    energy_contour += "=-x-"
            d["energy_contour"] = energy_contour

            mag_contour = "-o-"
            if len(image_list)==0:
                return None
            for i in xrange(1,len(image_list)):
                if(list_image_mags[i]>list_image_mags[i-1]):
                    mag_contour += "/-o-"
                elif list_image_mags[i]<list_image_mags[i-1]:
                    mag_contour += "\\-o-"
                else:
                    mag_contour += "=-o-"
            d["mag_contour"] = mag_contour

            # Difference between the first and maximum energies and magnetic moments
            deltaE_firstmax = max(list_image_energies) - list_image_energies[0]
            d["deltaE_firstmax"] = deltaE_firstmax
            deltaM_firstmax = max(list_image_mags) - list_image_mags[0]
            d["deltaM_firstmax"] = deltaM_firstmax

            # Difference between the last and maximum energies and magnetic moments
            deltaE_lastmax = max(list_image_energies) - list_image_energies[-1]
            d["deltaE_lastmax"] = deltaE_lastmax
            deltaM_lastmax = max(list_image_mags) - list_image_mags[-1]
            d["deltaM_lastmax"] = deltaM_lastmax

            # Difference between the endpoint energies and magnetic moments
            deltaE_endpoints = list_image_energies[-1] - list_image_energies[0]
            d["deltaE_endpoints"] = deltaE_endpoints
            deltaM_endpoints = list_image_mags[-1] - list_image_mags[0]
            d["deltaM_endpoints"] = deltaM_endpoints

            # Difference between the minimum and maximum energies and magnetic moments
            deltaE_maxmin = max(list_image_energies) - min(list_image_energies)
            d["deltaE_maxmin"] = deltaE_maxmin
            deltaM_maxmin = max(list_image_mags) - min(list_image_mags)
            d["deltaM_maxmin"] = deltaM_maxmin

            d["type"] = "NEB"

            return d
    """
