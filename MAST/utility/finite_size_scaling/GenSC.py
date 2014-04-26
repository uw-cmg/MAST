##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Wei Xie
# Last updated: 2014-04-25
##############################################################
import sys
import getopt 
import os
import shutil
import errno
#import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import pymatgen as mg
from pymatgen.analysis import ewald
from pymatgen.io import vaspio
import EneVsVm
from MAST.utility import MASTFile

def roof_mean(intA,intB):
    if (intA+intB)%2 == 0:
        return (intA+intB)/2
    else:
        return (intA+intB+1)/2

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
        
               

def minVertDist(inputlat):
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
               
def genLMN(primordial,minDefDist,maxNumAtoms):
    #LMN_UpLimit=int(maxNumAtoms//primordial.composition.num_atoms)
    
#    sc_list_raw=[]
    LMN_list_raw=[]
    Vm_list_raw=[]
    
    LMN_UpLimit=4
    for L in range(1,LMN_UpLimit+1):
            for M in range(1,LMN_UpLimit+1):
                    for N in range(1,LMN_UpLimit+1):
                        dummy_sc=primordial.copy()
                        dummy_sc.make_supercell([L,M,N])
                        if ((minVertDist(dummy_sc.lattice) >= minDefDist) and 
                         (dummy_sc.composition.num_atoms <= maxNumAtoms)):                                                       
                            #sc_list_raw.append(mg.core.structure.IStructure.from_sites(dummy_sc))
                            #sc_list_raw.append(dummy_sc)
                            LMN_list_raw.append([L,M,N])
                            Vm_list_raw.append(EneVsVm.CalcV_M(dummy_sc))

                            
    if len(LMN_list_raw) <= 5:
        LMN_list=LMN_list_raw
    else:
        LMN_list=[]
        Vm_LMN_dict={}
        for i in range(len(LMN_list_raw)):
            #print Vm_list_raw[i]
            #print LMN_list_raw[i]
            Vm_LMN_dict[Vm_list_raw[i]]=LMN_list_raw[i]
#            Vm_LMN_dict[str(EneVsVm.CalcV_M(sc_list_raw[i]))]=sc_list_raw[i]
#            Vm_LMN_dict[EneVsVm.CalcV_M(sc_list_raw[i])]=sc_list_raw[i]

        Vm_list_raw_sorted=sorted(Vm_list_raw)
        LMN_list.append(Vm_LMN_dict[Vm_list_raw_sorted[0]])
        LMN_list.append(Vm_LMN_dict[Vm_list_raw_sorted[
                               roof_mean(0,roof_mean(0,len(Vm_list_raw_sorted)-1))]])
        LMN_list.append(Vm_LMN_dict[Vm_list_raw_sorted[roof_mean(0,len(Vm_list_raw_sorted))]])
        LMN_list.append(Vm_LMN_dict[Vm_list_raw_sorted[roof_mean(len(Vm_list_raw_sorted),
                                                 roof_mean(0,len(Vm_list_raw_sorted)-1))]])
        LMN_list.append(Vm_LMN_dict[Vm_list_raw_sorted[len(Vm_list_raw_sorted)-1]])

    return LMN_list
   
        
def genSCs(primordial,lmn_list):
    cwdir = os.getcwd()
    genSCoutput = MASTFile()
    genSCoutput.data.append(cwdir+"\n")
    genSCoutput.data.append("----------------------------------\n")
    genSCoutput.data.append("ScalingLMN "+"    V_M     "+"Kpoint mesh" + '\n')
    genSCoutput.data.append("----------------------------------\n")
    print ("---------------------------------")
    print ("ScalingLMN " "    V_M     " "Kpoint mesh")
    print ("---------------------------------")
    for j in range(len(lmn_list)):
        dirname=str(lmn_list[j][0])+"x"+str(lmn_list[j][1])+"x"+str(lmn_list[j][2])
        #os.mkdir(dirname)
        mkdir_p(dirname)
        #shutil.copy("INCAR", dirname)
        #shutil.copy("POTCAR",dirname)
        os.chdir(dirname)
        dummystruct=primordial.copy()
        dummystruct.make_supercell(lmn_list[j])
        dummykpnt=mg.io.vaspio.Kpoints.automatic_density(dummystruct,1500)
        mg.write_structure(dummystruct,'POSCAR')
        dummykpnt.write_file('KPOINTS')
        print (str(lmn_list[j]) +" "+ str(EneVsVm.CalcV_M(dummystruct)) +" "+ str(dummykpnt.kpts[0]))
        genSCoutput.data.append((str(lmn_list[j]) +" "+ str(EneVsVm.CalcV_M(dummystruct)) +" "+ str(dummykpnt.kpts[0]))+"\n")
        os.chdir(cwdir)
    genSCoutput.to_file(os.path.join(cwdir,"supercell_list.txt"))
                
# def CalcNelect(defectCharge):  
"""
This is to calculate the NELECT parameter based ond defect charge and composition
"""

#Dirs = ["2x2x1","2x2x2","2x2x3","3x3x1","4x4x1","3x3x2"]
#for i in range(len(Dirs)):
#    PrLatMat=mg.io.vaspio.Poscar.from_file(Dirs[i]+"_VAl3+.vasp").structure.lattice
#    print Dirs[i]
#    print PrLatMat.a
#    print PrLatMat.b
#    print PrLatMat.c 
#    print minVertDist(PrLatMat)      

#testlat=mg.core.lattice.Lattice([[1,0,0],
#                                 [0.9,0.1,0],
#                                 [0.8,0.05,0.1]])
#print testlat.a #1.0
#print testlat.b #0.905538513814
#print testlat.c #0.80777472107
#print minVertDist(testlat) #0.141421356237


##kk=mg.io.vaspio.Poscar.from_file("POSCAR_pmcell_SiC").structure
#kk=mg.io.vaspio.Poscar.from_file("POSCAR").structure
#ll=genLMN(kk,4,1000)
##print ll
#genSCs(kk,ll)

