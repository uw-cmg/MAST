#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Wei Xie
# Last updated: 2014-05-14
##############################################################
import sys
import getopt 
import os
import shutil
import errno
import warnings
import numpy as np
from scipy import stats

import pymatgen as mg
from pymatgen.analysis import ewald
from pymatgen.io import vaspio

import EneVsVm 
import uniformPick    

class writer:
    """
    This is a helper class that copies all console ouput to a text file as well
    """
    def __init__(self, *writers) :
        self.writers = writers

    def write(self, text) :
        for w in self.writers:
            w.write(text)

def num_atoms_speciewise(inputstruct):
    """
    This is a helper method that output a list giving the number of atoms for 
    each atomic specie in the same order as POTCAR. 
    This method is used instead of pymatgen's composition class because the order 
    from pymatgen's composition class is not necessarify in the same order as POTCAR 
    """

    transitionSiteIndex=[0]
    SpecieLastSite=inputstruct.species[0]
    for i, specieThisSite in enumerate(inputstruct.species):
        if specieThisSite != SpecieLastSite:
            transitionSiteIndex.append(i)
        SpecieLastSite=specieThisSite
    numSpecieTypes=len(transitionSiteIndex)
    numAtomsSpeciewise=[None]*numSpecieTypes
    for j in range(numSpecieTypes-1):
        numAtomsSpeciewise[j]=transitionSiteIndex[j+1]-transitionSiteIndex[j]
    numAtomsSpeciewise[numSpecieTypes-1]=(int(inputstruct.composition.num_atoms)-
                                        transitionSiteIndex[numSpecieTypes-1])

    return numAtomsSpeciewise           
    
                        
def minVertDist(inputlat):
    """
    This is a helper method that calculates the the minimum distance between the
    vertices of an input 3-D lattice. Is is essentially the shortest distance 
    possible between two points in the space asuuming 3D periodic boundary conditions.
    Used to calculate the shortest defect-defect
    distance in a defected supercell containing one point defect. 
    """
    dummyMat = inputlat.matrix
    return min(np.linalg.norm(dummyMat[0]),
               np.linalg.norm(dummyMat[1]),
               np.linalg.norm(dummyMat[2]),
               np.linalg.norm(dummyMat[0]+dummyMat[1]),
               np.linalg.norm(dummyMat[0]+dummyMat[2]),
               np.linalg.norm(dummyMat[1]+dummyMat[2]),
               np.linalg.norm(dummyMat[0]-dummyMat[1]),
               np.linalg.norm(dummyMat[0]-dummyMat[2]),
               np.linalg.norm(dummyMat[1]-dummyMat[2]),
               np.linalg.norm(dummyMat[0]+dummyMat[1]+dummyMat[2]),
               np.linalg.norm(dummyMat[0]+dummyMat[1]-dummyMat[2]),
               np.linalg.norm(dummyMat[0]-dummyMat[1]+dummyMat[2]),
               np.linalg.norm(dummyMat[0]-dummyMat[1]-dummyMat[2])
               )
def defChg(inputStruct,inputPotcar,inputIncar):
    """
    This is a helper method that calculates the defect charge based on input files
    POSCAR, INCAR and POTCAR (not the file name, but objects of pymatgen's 
    corresponding Poscar, Incar, and Potcar class.
    """
    if ('NELECT' not in inputIncar): 
        raise RuntimeError("cannot find NELECT in input INCAR")
    natoms_el=num_atoms_speciewise(inputStruct)
    zval_el=[None]*len(inputPotcar)
    defchg=0
    for i in range(len(inputPotcar)):
        zval_el[i]=inputPotcar[i].ZVAL
        defchg=defchg+zval_el[i]*natoms_el[i]
    defchg=defchg-inputIncar['NELECT']
    return defchg  

def roof_mean(intA,intB):
    """
    This is a helper method that calculates the integer mean of two integers 
    rounded to roof if the direct mean of the two integers are not a integer.
    """
    if (intA+intB)%2 == 0:
        return (intA+intB)/2
    else:
        return (intA+intB+1)/2

#def genLMNs(primordial_struct,minDefDist=5,maxNumAtoms=600,minNumAtoms=64,numStructAsked=5):
def genLMNs(primordial_struct,minDefDist=8,maxNumAtoms=300,numStructAsked=6):
    """
    Calculates a list of LxMxN scaling factors for a input structure given inputs 
    of minmum defect-defect distance, maximum number of atoms that the scaled 
    supercell may have, and how many LxMxN scaling factors to outpout
    """   
    maxLMN=int((maxNumAtoms/primordial_struct.composition.num_atoms)**(1/3.0))
    print ("Generating candidate supercells with the scaling factors"+
           " L/M/N less than or equal to "+str(maxLMN)+".")
    print ("If you want to explore larger scaling factors L/M/N, "+ 
            "please rerun with larger maxNumAtoms input parameter.")
    print (" ")
    print maxLMN
    #####The following section generates a list of candidates LMN's######   
    Vm_LMN_dict={}
    print ("The following are candiate supercells:")   
    print ("ScalingLMN  " " V_M")
    print ("-----------------")
    for L in range(1,maxLMN+1):
            for M in range(1,maxLMN+1):
                    for N in range(1,maxLMN+1):
                        dummy_sc=primordial_struct.copy()
                        dummy_sc.make_supercell([L,M,N])
                        if ((minVertDist(dummy_sc.lattice) >= minDefDist) and 
                            #(dummy_sc.composition.num_atoms >= minNumAtoms) and
                            (dummy_sc.composition.num_atoms <= maxNumAtoms)):
                            dummy_Vm=np.asscalar(EneVsVm.CalcV_M(dummy_sc))

                            alreadyExist=False
                            for vmraw in Vm_LMN_dict.keys():
                                if abs(vmraw-dummy_Vm)<0.01:
                                    alreadyExist=True
                                    break                                                                                     
                            if not alreadyExist:
                                Vm_LMN_dict[dummy_Vm]=[L,M,N]
                                #print(str([L,M,N]) +"   "+ str(round(dummy_Vm,2)))
                                print(str([L,M,N]) +"   %4.2f" % dummy_Vm)
    print""   
    
    #####The following section selects "Good" LMN's from the candidates list######
    LMN_list_raw=Vm_LMN_dict.values()
    Vm_list_raw=Vm_LMN_dict.keys()
    
    if len(LMN_list_raw) < numStructAsked:
        LMN_list=LMN_list_raw
        warnings.warn("You asked for "+str(numStructAsked)+" supercells, but totally only "+str(len(LMN_list_raw))+" supercells were found to fulfill the input requirements!")
    elif len(LMN_list_raw) == numStructAsked:
        LMN_list=LMN_list_raw
    elif type(numStructAsked) != int or numStructAsked <=0:
        raise RuntimeError("You asked for "+str(numStructAsked)+ " supercells, \
            but this input parameter should be a positive integer number.")
    elif numStructAsked == 1:
        VmSelected=sorted(Vm_list_raw)[roof_mean(0,(len(Vm_list_raw)-1))]
        LMN_list=[Vm_LMN_dict[VmSelected]]
    elif numStructAsked == 2:
        VmSelected=[min(Vm_list_raw),max(Vm_list_raw)]
        LMN_list=[Vm_LMN_dict[VmSelected[0]],Vm_LMN_dict[VmSelected[1]]]                             
    else:
        VmSelected=uniformPick.pickSubList(Vm_list_raw,numStructAsked)        
        LMN_list=[]                
        for Vm in VmSelected:
            LMN_list.append(Vm_LMN_dict[Vm])

    print ("The following are selected supercells:")   
    print ("ScalingLMN  " " V_M")
    print ("-----------------")           
    for j in range(len(LMN_list)):
        dummystruct=primordial_struct.copy()
        dummystruct.make_supercell(LMN_list[j])
        print (str(LMN_list[j])+"   %4.2f" % EneVsVm.CalcV_M(dummystruct)) 
    print ""
    return LMN_list

def gensc(LMN_list,perf_primordial_struct, 
          def_primordial_struct,def_primordial_kpnt,
          def_primordial_potcar,def_primordial_incar):
    """
    Generate supercells given input LxMxN scaling factor list, the defected 
    primordial cell's POSCAR, INCAR, POTCAR and KPOINTS and undefected primordial
    cell's CONTCAR's corresponding object generated by pymatgen.
    """

    if ('MAGMOM' in def_primordial_incar):
        print ("WARNING: found MAGMOM in the primordial cell INCAR" 
                "please manually reset MAGMOM in the supercell INCAR "
                "to match it to the number of atoms. Otherwise the job will crash.")
    
    print""               
       
    cwdir = os.getcwd()
    print("To run the following supercells:")
#    print ("ScalingLMN " "  V_M    " " Kpoint  " " Natoms")
    print ("ScalingLMN " " Kpoint  "   "   Label")
    print ("---------------------------------") 
    for j in range(len(LMN_list)):
        dummystruct=perf_primordial_struct.copy()
        dummystruct.make_supercell(LMN_list[j])
        kppra=round((def_primordial_kpnt.kpts[0][0]*
                def_primordial_kpnt.kpts[0][1]*
                def_primordial_kpnt.kpts[0][2])/
                (1/def_primordial_struct.composition.num_atoms))        
        dummykpnt=mg.io.vaspio.Kpoints.automatic_density(dummystruct,kppra)
        dummykpnt.style=def_primordial_kpnt.style
        print (str(LMN_list[j])+"   " \
               #    +"%4.2f"      %EneVsVm.CalcV_M(dummystruct)+"   "+
              +str(dummykpnt.kpts[0][0])+"x"+ \
               str(dummykpnt.kpts[0][1])+"x"+ \
               str(dummykpnt.kpts[0][2])+"   "
               #+str(int(dummystruct.composition.num_atoms))
               +"label="
              +str(LMN_list[j][0])+"x"+ \
               str(LMN_list[j][1])+"x"+ \
               str(LMN_list[j][2])               
               )
 
        os.chdir(cwdir)
        
if __name__ == "__main__":
    saved = sys.stdout
    fout = file('out.log', 'w')
    sys.stdout = writer(sys.stdout, fout)
    
    perfDir = "perfect"
    defDir = "defect"
    
    ###read CONTCAR from primordial_perfect directory########
    cwdir=os.getcwd()
    os.chdir(perfDir)      
    perfectContcar=mg.io.vaspio.Poscar.from_file("CONTCAR").structure
    os.chdir(cwdir)

    #read CONTCAR, KPOINTS, POTCAR and INCAR from primordial_defect directory##
    os.chdir(defDir)       
    defectedContcar=mg.io.vaspio.Poscar.from_file("CONTCAR").structure
    defectedKpoints=mg.io.vaspio.Kpoints.from_file('KPOINTS')
    defectedPotcar=mg.io.vaspio.Potcar.from_file('POTCAR')
    defectedIncar=mg.io.vaspio.Incar.from_file('INCAR')
    os.chdir(cwdir)

    ##calculate and select a list of LMN scaling factors to run
    #os.chdir('..')
    LMN_list=genLMNs(perfectContcar)

    ##generate the input files for the above scaling factors
    gensc(LMN_list,perfectContcar,
          defectedContcar,defectedKpoints,defectedPotcar,defectedIncar)

    sys.stdout = saved
    fout.close()
