import os
import pymatgen
import numpy as np
from MAST.utility import MASTError


kboltz = 8.6173325E-5

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
            self.entrodict <dict of float>: freq-labeled dict of entropies
            self.formdict <dict of float>: freq-labeled dict of formation
            self.tempdict <dict of float>: temperature list of diffusion coeffs
    """
    def __init__(self, directory, tempstart=73, tempend=1273, tempstep=100, freqmodel=0, freqdict=None):
        """
            Args:
                directory <str>: directory to start in
                tempstart <float>: Temp in K to start range (default 73)
                tempend <float>: Temp in K to end range (default 1273)
                tempstep <float>: step in K for temp. range (default 100)
                freqmodel <int>: Integer for frequency model
                freqdict <dict>: Dictionary for frequency-to-hop-labels
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
        try:
            pm = PickleManager(self.directory + 'input_options.pickle')
            self.input_options = pm.load_variable()
        except OSError:
            pass
        return
        self.set_up_temp_dict(tempstart, tempend, tempstep)
        self.set_freqmodel(freqmodel)
        self.set_frequency_dict(freqdict)

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
        return self.tempdict

    def set_freqmodel(self, freqmodel=0)        
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
            raise MASTError(self.__class__.__name__, "%s is not a supported frequency model." % trymodel)
        self.freqmodel = trymodel
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
        return self.freqdict


    def set_host_and_solute(self):
        """Find the host and the solute based on relative numbers in POSCAR,
            for a BINARY, DILUTE system. (One host element type, 
            one solute element type.)
        """
        label=""
        if self.freqmodel == 1:
            self.solute=None
            if not 'w0' in self.freqdict.keys():
                raise MASTError(self.__class__.__name__, "Need a label for frequency w0 for one-frequency model.")
            label = self.freqdict['w0']
        else:
            if not 'w2' in self.freqdict.keys():
                raise MASTError(self.__class__.__name__, "Not enough frequencies to determine solute and host.")
            label = self.freqdict['w2']       
        [startlabel, tstlabel]=self.find_start_and_tst_labels(label)
        mydir = self.get_labeled_dir(startlabel)
        pospath = os.path.join(mydir, "POSCAR")
        if not os.path.isfile(pospath):
            raise MASTError(self.__class__.__name__, "No POSCAR in %s." % pospath)
        mypos = pymatgen.io.vaspio.Poscar.from_file(pospath)
        sitesym = mypos.site_symbols
        natoms = mypos.natoms
        hostidx=natoms.index(max(natoms))
        self.host=sitesym[hostidx]
        if self.freqmodel > 1:
            solidx=natoms.index(min(natoms))
            self.solute=sitesym[solidx]
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
            if (label in myname) and (app1 in myname) and (app2 in myname):
                tstlabel = label
                startlabel = label.split('-')[0]
                break
            elif (backlab in myname) and (app1 in myname) and (app2 in myname):
                tstlabel = label
                startlabel = label.split('-')[1]
                break
        if tstlabel == "":
            raise MASTError(self.__class__.__name__, 
                "No directory found for hop label %s." % label)
        return [startlabel, tstlabel]

    def get_labeled_dir(self, label, append="stat"):
        """Get a labeled directory.
            Args:
                label <str>: label for the directory, like "vac1"
                append <str>: extra piece which should appear in the 
                                directory name
            Returns:
                mypath <str>: The FIRST matching directory found, from a 
                                sorted list. Returns a full path.
        """
        dirlist = os.listdir(self.directory)
        dirlist.sort()
        mydir=""
        for mydir in dirlist:
            if (label in mydir) and (append in mydir):
                mypath=os.path.join(self.directory, mydir)
                break
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
                <float>: transition state energy - starting state energy
        """
        [startlabel, tstlabel] = find_start_and_tst_labels(stem, label)
        startenergy = self.get_labeled_energy(stem, startlabel)
        maxenergy = self.get_max_energy(stem, tstlabel)
        if (startenergy > maxenergy):
            print "Attention! Maximum energy is from an endpoint, not an image!"
        return (maxenergy - startenergy)

    def get_labeled_energy(self, label, append="stat"):
        """Get an energy from a directory.
            Args:
                label <str>: label for directory
                append <str>: additional tag for the directory
        """
        mypath = self.get_labeled_dir(label, append)
        myenergy=0
        myenergy = self.get_total_energy(mypath)
        return myenergy

    def get_total_energy(directory): #synchronize this call with Glen's 
        """Returns the total energy from a directory.
            Args:
                directory <str>: directory (full path)
            Returns:
                <float>: Total energy (E0 from VASP)
            SHOULD BE REPLACED WITH PROGRAM-INDEPENDENT SWITCH.
        """
        from pymatgen.io.vaspio.vasp_output import Vasprun
        # Are we dealing with VASP?
        if ('vasprun.xml' in os.listdir(directory)):
        # Modified from the PyMatGen Vasprun.final_energy() function to return E_0
            return Vasprun('%s/vasprun.xml' % directory).ionic_steps[-1]["electronic_steps"][-1]["e_0_energy"]

    def get_max_energy(self, tstlabel):
        """Get maximum image energy from an NEB directory.
            Args:
                tstlabel <str>: label for the NEB directory
            Returns:
                <float>: Maximum energy among the images.
        """
        mypath = self.get_labeled_dir(tstlabel)
        imdirs = os.listdir(mypath)
        imdirs.sort()
        mydir=""
        imenergs=list()
        for mydir in imdirs:
            imenergs[self.get_total_energy(os.path.join(mypath,mydir))] = mydir
        maxenerg = max(imenergs)
        maxdir = imenergs[maxenerg]
        return maxenerg
   

    def get_hopdict(self):
        """Get dictionary of hops."""
        freqs = self.freqdict.keys()
        hopdict=dict()
        for freq in freqs:
            neblabel = freqdict[freq]
            hopdict[freq] = self.neb_barrier(neblabel)
        self.hopdict = hopdict
        return self.hopdict

    def get_formdict(self):
        """Get formation energy dictionary.
        """
        formdict=dict()
        #pure_def = #TTM DEBUG NEED TO SWITCH TO USE GLEN'S FUNCTION
        for freq in self.freqdict.keys():
            formdict[freq] = self.get_defect_energy_simple(freqdict[freq])
        self.formdict = formdict
        return self.formdict

    def get_defect_energy_simple(self, label, append="stat", perf="perf"):
        """Need to replace with Glen's function.
            Get simplified defect formation energy. NO CHARGED CELLS.
            Args:
                label <str>: label for hop
                append <str>: append for name 
                perf <str>: perfect-cell designation
            Return:
                <float>: defect formation energy
        """
        [startlabel, tstlabel] = self.find_start_and_tst_labels(label)
        defpath=self.get_labeled_dir(startlabel, append)
        perfpath=self.get_labeled_dir(perf, append)
        defectedenergy = self.get_total_energy(defpath)
        perfectenergy = self.get_total_energy(perfpath)
        if defectedenergy==None or perfectenergy==None:
            raise MASTError(self.__class__.__name__,"No defected energy found.")
        delta = defectedenergy - perfectenergy
        ecohpure=0 # Need to get this from a calc??
        if ecohpure==0:
            mypos=pymatgen.io.vaspio.Poscar.from_file(os.path.join(perfpath, "POSCAR"))
            totatoms=sum(mypos.natoms)
            echopure=perfectenergy/totatoms
            myfile=MASTFile()
            myfile.data.append("WARNING: Cohesive energy has not been properly calculated in %s because a pseudopotential reference has not been taken into account." % dirname)
            myfile.to_file(os.path.join(dirname,"WARNING"))
        eform = delta + ecohpure #Add on reference atom
        return eform

    def get_attemptdict(self):
        """Get approximate attempt vibrational frequency dictionary, 
            using Adams approximations:
            v_0=v_1=v_3=v_4 and (v2/v0)=sqrt(mass0*Tmeltpure2/(mass2*Tmeltpure0)
            Temperatures are in KELVIN.
            Masses are in AMU.
        """
        stock_v=1e13
        vibdict=dict()
        hostelem = pymatgen.core.periodic_table.Element(self.host)
        solelem = pymatgen.core.periodic_table.Element(self.solute)
        masshost = hostelem.atomic_mass
        masssolute = solelem.atomic_mass
        tmeltpurehost = self.parse_melting_point(hostelem)
        tmeltpuresolute = self.parse_melting_point(solelem)
        freqs = freqdict.keys()
        freqs.sort()
        for freq in freqs:
            if not (freq == 'w2'):
                vibdict[freq] = stock_v
            else:
                vibdict[freq] = vibdict['w0'] * np.sqrt(masshost*tmeltpuresolute/(masssolute*tmeltpurehost))
        self.attemptdict=vibdict
        return self.attemptdict
    

    def get_entrodict_approx(self):
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
        entrodict=dict()
        stock_S_v=2.07e-4 #Shewmon, Sv/R approx 2.4?? =2.07e-4
        stock_S_m=2.07e-4
        beta=0.4
        hostelem = pymatgen.core.periodic_table.Element(self.host)
        tmeltpurehost = self.parse_melting_point(hostelem)
        for freq in freqdict.keys():
            if not (freq == 'w2'):
                entrodict[freq] = stock_S_v + stock_S_m
            else:
                eact2 = self.hopdict['w2'] + self.formdict['w2']
                eact1 = self.hopdict['w1'] + self.formdict['w1']
                delta_act_s = beta*(eact2-eact1)/tmeltpurehost
                entrodict['w2'] = entrodict['w0'] + delta_act_s
        self.entrodict = entrodict
        return self.entrodict

    def parse_melting_point(self, elem):
        """Parse the melting point.
            Args:
                elem <Element>: pymatgen Element class object
            Return:
                <float>: melting point in Kelvin
        """
        meltstr = str(elem.melting_point) #str(<unicode>)
        meltval = meltstr.split()[0]
        return float(meltval)


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
        self.get_entrodict_approx()
        self.get_vibdict_approx() #Adams approximations
        freqs = self.freqdict.keys()
        jumpfreqdict=dict()
        for freq in freqs:
            jumpfreqdict[freq] = self.estimate_jump_freq(self.vibdict[freq],
                self.entrodict[freq],
                self.hopdict[freq]+self.formdict[freq], temp)
        jfw0 = jumpfreqdict['w0']
        jfw1 = jumpfreqdict['w1']
        jfw2 = jumpfreqdict['w2']
        jfw3 = jumpfreqdict['w3']
        jfw4 = jumpfreqdict['w4']
        f_0 = 0.715
        myx = jfw4/jfw0
        bigf_num = 10*np.power(myx,4) + 180.5*np.power(myx,3) + 927*np.power(myx,2) + 1341
        bigf_denom = 2*np.power(myx,4) + 40.2*np.power(myx,3) + 254*np.power(myx,2) + 597*myx + 435
        bigf_x = 1-(1/7)*(bigf_num/bigf_denom)
        f_2_num = 1+3.5*bigf_x *(jfw3/jfw1)
        f_2_denom = 1+(jfw2/jfw1) + 3.5*bigf_x*(jfw3/jfw1)
        f_2 = f2_num / f_2_denom
        vacconc = self.get_vac_conc()
        lattparam = self.get_latt_param()
        Dself = vacconc*lattparam^2*jfw0
        Dsol_num = Dself * f_2 * jfw2 * jfw4 * jfw1
        Dsol_denom = f_0 * jfw1 * jfw0 * jfw3
        Dsolute = Dsol_num / Dsol_denom
        return [Dself, Dsolute]

    def get_vac_conc(self, perf="perf"):
        """Get vacancy concentration.
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
        return vacconc


    def get_latt_param(self, perf="perf"):
        """Get lattice parameter
            Args:
                perf <str>: designation string for perfect cell
                            (default) "perf"
            Returns:
                <float>: lattice parameter in Angstroms.
        """
        perfpath=self.get_labeled_dir(perf)
        perfpospath = os.path.join(perfpath, "POSCAR")
        if not os.path.isfile(perfpospath):
            raise MASTError(self.__class__.__name__, "No POSCAR in %s." % perfpospath)
        perfpos = pymatgen.io.vaspio.Poscar.from_file(perfpospath)
        perfstr = perfpos.structure
        mylatt = perfstr.lattice.a
        return mylatt
    
    def estimate_jump_freq(self, vibfreq, actentropy, actenergy, temp):
        """Estimate jump frequency (Adams, Foiles, Wolfer)
            w_i = v_i exp(S_i/k) exp(-E_i/kT)
        """
        jumpfreq = vibfreq*np.exp(actentropy/kboltz)*np.exp(-1*actenergy/(kboltz*temp))
        return jumpfreq

    def diffusion_coefficient(self):
        """Diffusion coefficient switcher and temperature looper."""
        tkeys = self.tempdict.keys()
        for tkey in tkeys:
            if self.freqmodel == 1:
                self.tempdict[tkey] = self.one_freq(tkey)
            elif self.freqmodel == 5:
                self.tempdict[tkey] = self.five_freq(tkey)
            else:
                raise MASTError(self.__class__.__name__, "%s is not a supported frequency model." % self.freqmodel)
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
        self.get_entrodict_approx()
        self.get_vibdict_approx() #Adams approximations
        freqs = self.freqdict.keys()
        jumpfreqdict=dict()
        for freq in freqs:
            jumpfreqdict[freq] = self.estimate_jump_freq(self.vibdict[freq],
                self.entrodict[freq],
                self.hopdict[freq]+self.formdict[freq], temp)
        jfw0 = jumpfreqdict['w0']
        f_0 = 0.715
        vacconc = self.get_vac_conc()
        lattparam = self.get_latt_param()
        Dself = vacconc*lattparam^2*jfw0
        return Dself

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
    directory=os.path.expandvars(sys.argv[1])
    tstart=73
    tend=1273
    tstep=100
    freqmodel=1
    try:
        trytstart=sys.argv[2]
        tstart=trytstart
    except IndexError:
        pass
    try:
        trytend=sys.argv[3]
        tend=trytend
    except IndexError:
        pass
    try:
        trytstep=sys.argv[4]
        tstep=trytstep
    except IndexError:
        pass
    try:
        tryfreqmodel=sys.argv[5]
        freqmodel=tryfreqmodel
    except IndexError:
        pass
    freqdict=None
    try:
        fdstr = sys.argv[6]
        freqdict=dict()
        fdlist=fdstr.strip()split(",")
        for fditem in fdlist:
            mykey = fditem.split("=")[0]
            myval = fditem.split("=")[1]
            freqdict[mykey]=myval
    except IndexError:
        pass

    print 'Looking at diffusion coefficient for %s' % directory
    DC = DiffusionCoefficient(directory, tstart, tend, tstep, freqmodel, freqdict)
    print DC.diffusion_coefficient()

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
