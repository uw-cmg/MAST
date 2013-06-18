import os
import pymatgen
import numpy as np
from MAST.utility import MASTError


kboltz = 8.6173325E-5

def diffusion_coefficient(stem="", freqmodel=5, freqdict=dict(), temp=1173, vacconc, lattparam):
    """Wrapper to dilute solute diffusion coefficient calculations.
        DOES NOT SUPPORT CHARGED SUPERCELLS YET.
        Args:
            stem <str>: stem of system/recipe to calculate (relative to mydir)
            freqmodel <int>: frequency model   
                                5 (default)
                                1 (one frequency; pure system)
            freqdict <dict of str>: frequency-to-label dictionary:
                        freqdict['w0']=['vac1-vac2']
                        freqdict['w4']=['vac10-vac9']
                            If a label is backwards from what is in the 
                            neb ingredient's actual name, then
                            the backwards hop will be taken.
                            Initial state for phonons will be the first
                            label in the label pair.
            temp <int>: temperature in Kelvin
                            1173 (K, default)
        Returns:
            D_0 and Q, and D at the temperature specified
            for D* = D_0 * exp(-Q/kT)
    """
    if stem == "":
        raise MASTError("postprocessing, diffusion_coefficient",
            "Recipe stem was not given. Erroring out.")
        
    if freqmodel == 5:
        [myD0,myQ]=five_freq(stem, freqdict)
    elif freqmodel == 1:
        [myD0,myQ]=one_freq(stem, freqdict)
    else:
        raise MASTError("postprocessing, diffusion_coefficient",
            "Frequency model %s is not supported." % freqmodel)
    myD = myD0 * np.exp(-1*myQ/(kboltz*temp))
    return myD

def neb_barrier(stem, label):
    """Get the NEB barrier."""
    [startlabel, tstlabel] = find_start_and_tst_labels(stem, label)
    startenergy = get_labeled_energy(stem, startlabel)
    maxenergy = get_max_energy(stem, tstlabel)
    if (startenergy > maxenergy):
        print "Attention! Maximum energy is from an endpoint, not an image!"
    return (maxenergy - startenergy)

def get_labeled_energy(stem, label, append="stat"):
    """Get an energy from a directory."""
    dirname = os.path.dirname(stem)
    dirlist = os.path.listdir(dirname)
    dirlist.sort()
    mydir=""
    myenergy=0
    for mydir in dirname:
        if (label in mydir) and (append in mydir):
            myenergy = get_total_energy(mydir)
            break
    return myenergy

def get_max_energy(stem, tstlabel, imageflag="image"):
    """Get maximum image energy from an NEB directory."""
    dirname = os.path.dirname(stem)
    dirlist = os.path.listdir(dirname)
    dirlist.sort()
    mydir=""
    imdirs=list()
    for mydir in dirlist:
        if (imageflag in mydir) and (tstlabel) in mydir:
            imdirs.append(mydir)
    imenergs=dict()
    mydir=""
    for mydir in imdirs:
        imenergs[get_total_energy(mydir)] = mydir
    maxenerg = max(imenergs)
    maxdir = imenergs[maxenerg]
    return maxenerg
    

def get_total_energy(directory): #synchronize this call with Glen's 
    """Returns the total energy from a directory"""
    from pymatgen.io.vaspio.vasp_output import Vasprun
    # Are we dealing with VASP?
    if ('vasprun.xml' in os.listdir(directory)):
    # Modified from the PyMatGen Vasprun.final_energy() function to return E_0
        return Vasprun('%s/vasprun.xml' % directory).ionic_steps[-1]["electronic_steps"][-1]["e_0_energy"]

def find_start_and_tst_labels(stem, label):
    """Finds the start label and transition state labels for ingredients.
        Args:
            stem <str>: stem name
            label <str>: neb label 
        Returns:
            [startlabel, tstlabel]
            startlabel <str>: label for start directory (e.g. "vac1" or "vac2")
            tstlabel <str>: label for neb (e.g. "vac1-vac2" or "vac2-vac1")
    """
    startlabel=""
    tstlabel=""
    mydir = os.path.dirname(stem)
    mybase = os.path.basename(stem)
    backlabel = reverse_label(label)
    for myname in os.listdir(mydir):
        if (label in myname) and (mybase in myname):
            tstlabel = label
            startlabel = label.split('-')[0]
            break
        elif (backlabel in myname) and (mybase in myname):
            tstlabel = label
            startlabel = label.split('-')[1]
            break
    if tstlabel == "":
        raise MASTError("utility, diffusioncoefficient", 
            "No directory found for hop label %s." % label)
    return [startlabel, tstlabel]

def reverse_label(label):
    """Reverse a label"""
    lsplit = label.split('-')
    backlabel = '-'.join([lsplit[1],lsplit[0]])
    return backlabel

def neb_profile():
    """Get the NEB profile."""
    pass

def one_freq(stem, freqdict):
    """Pure diffusion, one frequency
        Args:
            (see diffusion_coefficient)
    """
    pass

def five_freq(stem, freqdict, temp, vacconc, lattparam, tmeltpurehost, tmeltpuresolute, masshost, masssolute, ecohpure, sys="FCC"):
    """Five-frequency model.
        Args:
            (see diffusion_coefficient)
    Adams, J. B., Foiles, S. M. & Wolfer, W. G. 
    Self-diffusion and impurity diffusion of fcc metals 
    using the five-frequency model and the Embedded Atom Method. 
    Journal of Materials Research 4, 102.112 (2011)
    """
    hopdict = get_hopdict(stem, freqdict)
    formdict = get_formdict(stem, freqdict, ecohpure)
    entrodict = get_entrodict_approx(stem, freqdict, hopdict, formdict, tmeltpurehost) #Adams approximations
    vibdict = get_vibdict_approx(stem, freqdict, tmeltpurehost, tmeltpuresolute, masshost, masssolute) #Adams approximations
    freqs = freqdict.keys()
    jumpfreqdict=dict()
    for freq in freqs:
        jumpfreqdict[freq] = estimate_jump_freq(vibdict[freq],entrodict[freq],
                hopdict[freq]+formdict[freq], temp)
    jfw0 = jumpfreqdict['w0']
    jfw1 = jumpfreqdict['w1']
    jfw2 = jumpfreqdict['w2']
    jfw3 = jumpfreqdict['w3']
    jfw4 = jumpfreqdict['w4']
    if not sys=="FCC":
        raise MASTError("utility, diffusioncoefficient", "Only calculates five-frequency model for FCC.")
    if sys=="FCC":
        print "FCC:"
        f_0 = 0.715
        myx = jfw4/jfw0
        bigf_num = 10*np.power(myx,4) + 180.5*np.power(myx,3) + 927*np.power(myx,2) + 1341
        bigf_denom = 2*np.power(myx,4) + 40.2*np.power(myx,3) + 254*np.power(myx,2) + 597*myx + 435
        bigf_x = 1-(1/7)*(bigf_num/bigf_denom)
        f_2_num = 1+3.5*bigf_x *(jfw3/jfw1)
        f_2_denom = 1+(jfw2/jfw1) + 3.5*bigf_x*(jfw3/jfw1)
        f_2 = f2_num / f_2_denom
        Dself = vacconc*lattparam^2*jfw0
        Dsol_num = Dself * f_2 * jfw2 * jfw4 * jfw1
        Dsol_denom = f_0 * jfw1 * jfw0 * jfw3
        Dsolute = Dsol_num / Dsol_denom
    return [Dself, Dsolute]

def estimate_jump_freq(vibfreq, actentropy, actenergy, temp):
    """Estimate jump frequency (Adams, Foiles, Wolfer)
        w_i = v_i exp(S_i/k) exp(-E_i/kT)
    """
    jumpfreq = vibfreq*np.exp(actentropy/kboltz)*np.exp(-1*actenergy/(kboltz*temp))
    return jumpfreq
    

def get_hopdict(stem, freqdict):
    """Get dictionary of hops."""
    freqs = freqdict.keys()
    hopdict=dict()
    for freq in freqs:
        neblabel = freqdict[freq]
        hopdict[freq] = neb_barrier(stem, neblabel)
    return hopdict

def get_formdict(stem, freqdict, ecohpure):
    """Get formation energy dictionary.
    """
    formdict=dict()
    pure_def = #TTM DEBUG NEED TO SWITCH TO USE GLEN'S FUNCTION
    for freq in freqdict.keys():
        formdict[freq] = get_defect_energy_simple(stem, freqdict[freq], ecohpure)
    return formdict

def get_defect_energy_simple(stem, label, append="stat", perf="perf", ecohpure):
    """Need to replace with Glen's function"""
    [startlabel, tstlabel] = find_start_and_tst_labels(stem, freqdict[freq])
    dirname = os.path.dirname(stem)
    dirlist = os.listdir(dirname)
    dirlist.sort()
    mydir=""
    defectedenergy=None
    perfectenergy=None
    for mydir in dirlist:
        if (startlabel in mydir) and (append in mydir):
            defectedenergy = get_total_energy(mydir)
        elif (perf in mydir) and (append in mydir):
            perfectenergy = get_total_energy(mydir) 
    if defectedenergy==None or perfectenergy==None:
        raise MASTError("utility diffusioncoefficient","No defected energy found.")
    delta = defectedenergy - perfectenergy
    eform = delta + ecohpure #Add on reference atom
    return eform

def get_entrodict_approx(stem, freqdict, hopdict, formdict, tmeltpurehost):
    """Get approximate entropy.
        Use Adams approximations for now, Wert-Zener, where S_0=S_1=S_3=S_4 and
        S_2 - S_0 = beta(E_2 - E_1)/Tmelt, pure,
        where S are activation entropies and
        E are activation energies and beta is assumed to be 0.4

        Dobson et al, S/k approx 10?            
        Shewmon, Sm approx beta*Q/Tm, Sv/R = 2.4
        http://www.tf.uni-kiel.de/matwis/amat/
            def_en/kap_3/backbone/r3_1_1.html, Sform ~ Smig ~ 1k?         
    """
    entroform_d=dict()
    entromig_d=dict()
    entrodict=dict()
    S_v=2.07e-4
    beta=0.4
    for freq in freqdict.keys():
        if not (freq == 'w2'):
            entroform_d[freq] = S_v #Shewmon, Sv/R approx 2.4?? =2.07e-4
            entromig_d[freq] = S_v
            entrodict[freq] = entroform_d[freq] + entromig_d[freq]
        else:
            eact2 = hopdict['w2'] + formdict['w2']
            eact1 = hopdict['w1'] + formdict['w1']
            delta_act_s = beta*(eact2-eact1)/tmeltpurehost
            entrodict['w2'] = entrodict['w0'] + delta_act_s
    return entrodict

def get_vibdict_approx(stem, freqdict, tmeltpurehost, tmeltpuresolute, masshost, masssolute):
    """Get approximate vibrational dictionary, using Adams approximations:
        v_0=v_1=v_3=v_4 and (v2/v0) = sqrt(mass0*Tmeltpure2/(mass2*Tmeltpure0)
    """
    stock_v=1e13
    vibdict=dict()
    for freq in freqdict.keys():
        if not (freq == 'w2'):
            vibdict[freq] = stock_v
        else:
            vibdict[freq] = vibdict['w0'] * np.sqrt(masshost*tmeltpuresolute/(masssolute*tmeltpurehost))
    return vibdict
    

def phonons():
    """Get phonon info???"""
    pass

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
