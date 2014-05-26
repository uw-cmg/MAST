#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Wei Xie
# Last updated: 2014-04-19
##############################################################
import sys
import getopt 
import os
import shutil
import errno
import numpy as np
from scipy import stats

import pymatgen as mg
from pymatgen.analysis import ewald
from pymatgen.io import vaspio

import EneVsVm     

class writer:
        def __init__(self, *writers) :
                self.writers = writers

        def write(self, text) :
                for w in self.writers :
                        w.write(text)

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: 
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise      
def num_atoms_speciewise(inputstruct):

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
    
def defChg(inputStruct,inputPotcar,inputIncar):
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

def roof_mean(intA,intB):
    if (intA+intB)%2 == 0:
        return (intA+intB)/2
    else:
        return (intA+intB+1)/2

               
def list2Bins(inputlist,numBins):
    '''
    helper function to group the elements of inputlist into some evenly spaced 
    subgroups called bins. Each bin has two boundaries and there are totally numBins+1
    bin boundaries.    
    '''   
    listMin=min(inputlist)
    listMax=max(inputlist)

    if (numBins <= 0):
        raise RuntimeError("numBins <= 0")
    elif (listMin == listMax):
        raise RuntimeError("All memebers of the list have the same value")
    else:
        binWidth=(listMax-listMin)/numBins

        binLowerBounds=np.arange(listMin,listMax,binWidth)
        #binUpperBounds=binLowerBounds+binWidth
        binBounds=np.append(binLowerBounds,listMax)
        
        binLowerBounds=binLowerBounds.tolist()
        #binUpperBounds=binUpperBounds.tolist()
        binBounds=binBounds.tolist()     

        #itemsClosest2binBounds=[[] for i in range(len(binBounds))]
        #itemsClosest2binBounds_index=[[] for i in range(len(binBounds))]
        #closestDist2binBounds=[[] for i in range(len(binBounds))]
 
        itemsClosest2binBounds=[None]*(len(binBounds))
        #itemsClosest2binBounds_index=[None]*(len(binBounds))
        closestDist2binBounds=[None]*(len(binBounds))
        
        for i, item in enumerate(inputlist):
            for j, boundary in enumerate(binBounds):
                dummyDist=item-boundary
                if (-binWidth/2<=dummyDist) and (dummyDist<binWidth/2):
                   if ((itemsClosest2binBounds[j]==None) or 
                       (abs(dummyDist)<closestDist2binBounds[j])):
                       itemsClosest2binBounds[j]=inputlist[i]
                       #itemsClosest2binBounds_index[j]=[i]
                       closestDist2binBounds[j]=abs(dummyDist)
                   break                       
                                              
    #return (binBounds,itemsClosest2binBounds,closestDist2binBounds)    
    #return (itemsClosest2binBounds,itemsClosest2binBounds_index)
    return itemsClosest2binBounds


#def genLMNs(primordial_struct,minDefDist=5,maxNumAtoms=600,minNumAtoms=64,numStructAsked=5):
def genLMNs(primordial_struct,minDefDist=5,maxNumAtoms=600,numStructAsked=5):   
    maxLMN=int((maxNumAtoms/primordial_struct.composition.num_atoms)**(1/3.0))
    print ("Generating candidate supercells with the scaling factors"+
           " L/M/N less than or equal to "+str(maxLMN)+".")
    print ("If you want to explore larger scaling factors L/M/N, " + 
           "please rerun with larger maxNumAtoms input parameter.")
    print (" ")
    print maxLMN
    #####The following section generate a list of candidates LMN's######   
    LMN_list_raw=[]
    Vm_list_raw=[]
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
                            #print str([L,M,N]),"tmp"                            
                            alreadyExist=False
                            for vmraw in Vm_list_raw:
                                #print str([L,M,N]),"good"
                                if abs(vmraw-dummy_Vm)<0.00001:
                                    alreadyExist=True
                                    break                                                                                     
                            #if dummy_Vm not in Vm_list_raw:
                            if not alreadyExist:
                                LMN_list_raw.append([L,M,N])
                                Vm_list_raw.append(dummy_Vm)
                                print str([L,M,N]),"start"
                                print(str([L,M,N]) +"   "+ str(round(dummy_Vm,2)))
                                print str([L,M,N]),"finish"  
    print""   
    
    #####The following section select "Good" LMN's from the candidates list######
    if len(LMN_list_raw) < numStructAsked:
        LMN_list=LMN_list_raw
        print ('Warning: You asked for '+str(numStructAsked)+ 
               ' supercells, but totally only '
               +str(len(LMN_list_raw))+
               ' supercells were found to fullfile the input requirements!')
    elif len(LMN_list_raw) == numStructAsked:
        LMN_list=LMN_list_raw
    else:
        LMN_list=[]

        Vm_LMN_dict={}
        for i in range(len(LMN_list_raw)):
            Vm_LMN_dict[Vm_list_raw[i]]=LMN_list_raw[i]      
        
        Vm_list_raw_copy=list(Vm_list_raw)
        
        VmClosest2binBounds=list2Bins(Vm_list_raw,numStructAsked-1)
        
        VmSelected=[]        
        for kk in VmClosest2binBounds:
            if kk!=None:
                VmSelected.append(kk)                                              
                Vm_list_raw_copy.remove(kk)
                
        while len(VmSelected) < numStructAsked:                       
            numStructDeficit=numStructAsked-len(VmSelected)

            if (numStructDeficit == 1):
                Vm_list_raw_copy_sorted=sorted(Vm_list_raw_copy)
                dummyVm=Vm_list_raw_copy_sorted[roof_mean(0,len(Vm_list_raw_copy_sorted))]
                VmSelected.append(dummyVm)                                              
                #Vm_list_raw_copy.remove(dummyVm)
            else:
                VmClosest2binBounds_loop=list2Bins(Vm_list_raw_copy,numStructDeficit-1)
                for mm in VmClosest2binBounds_loop:
                    if mm!=None:
                        VmSelected.append(mm)
                        Vm_list_raw_copy.remove(mm)
                        
        for Vm in sorted(VmSelected):
            LMN_list.append(Vm_LMN_dict[Vm])

    print ("The following are selected supercells:")   
    print ("ScalingLMN  " " V_M")
    print ("-----------------")           
    for j in range(len(LMN_list)):
        dummystruct=primordial_struct.copy()
        dummystruct.make_supercell(LMN_list[j])
        print (str(LMN_list[j])+"   "+
               str(round(EneVsVm.CalcV_M(dummystruct),2))) 
    print ""
    return LMN_list

def gensc(LMN_list,perf_primordial_struct, 
          def_primordial_struct,def_primordial_kpnt,
          def_primordial_potcar,def_primordial_incar):

    if ('MAGMOM' in def_primordial_incar):
        print ("WARNING: found MAGMOM in the primordial cell INCAR" 
                "please manually reset MAGMOM in the supercell INCAR "
                "to match it to the number of atoms. Otherwise the job will crash.")
    
    defchg=defChg(def_primordial_struct,def_primordial_potcar,def_primordial_incar)
    print("The defect charge is "+str(defchg)+" in the defected primordial cell.")
    print""               
       
    cwdir = os.getcwd()
    print("To run the following supercells:")
    print ("ScalingLMN " "  V_M    " " Kpoint  " " Natoms")
    print ("-------------------------------------") 
    for j in range(len(LMN_list)):

        dirname=str(LMN_list[j][0])+"x"+str(LMN_list[j][1])+"x"+str(LMN_list[j][2])
        mkdir_p(dirname)
        os.chdir(dirname)
        
        dummystruct=perf_primordial_struct.copy()
        dummystruct.make_supercell(LMN_list[j])
        ####The next section does embeding#######            

        ####This above section does embeding#######
        mg.write_structure(dummystruct,'POSCAR') 
       
        kppra=round((def_primordial_kpnt.kpts[0][0]*
                def_primordial_kpnt.kpts[0][1]*
                def_primordial_kpnt.kpts[0][2])/
                (1/def_primordial_struct.composition.num_atoms))        
        dummykpnt=mg.io.vaspio.Kpoints.automatic_density(dummystruct,kppra)
        dummykpnt.style=def_primordial_kpnt.style
        dummykpnt.write_file('KPOINTS')
                
        dummystruct_natoms_el=num_atoms_speciewise(dummystruct)
        dummynelect=-defchg
        for i in range(len(def_primordial_potcar)):
            dummyZval=def_primordial_potcar[i].ZVAL
            dummynelect=dummynelect+dummyZval*dummystruct_natoms_el[i]
        def_primordial_incar.__setitem__('NELECT',dummynelect)
        def_primordial_incar.write_file('INCAR')
        
        def_primordial_potcar.write_file('POTCAR')
        print (str(LMN_list[j])+"   "+
               str(round(EneVsVm.CalcV_M(dummystruct),2))+"   "+
               str(dummykpnt.kpts[0])+"   "+
               str(int(dummystruct.composition.num_atoms)))
 
        os.chdir(cwdir)
        
if __name__ == "__main__":
    saved = sys.stdout
    fout = file('out.log', 'w')
    sys.stdout = writer(sys.stdout, fout)
    
    perfDir = "primordial_perfect"
    defDir = "primordial_defected"
    

    jj=mg.io.vaspio.Poscar.from_file(perfDir+"_CONTCAR").structure 

    kk=mg.io.vaspio.Poscar.from_file(defDir+"_CONTCAR").structure
    ll=mg.io.vaspio.Kpoints.from_file(defDir+"_KPOINTS")
    mm=mg.io.vaspio.Potcar.from_file(defDir+'POTCAR')
    nn=mg.io.vaspio.Incar.from_file(defDir+'INCAR')

    LMN_list=genLMNs(jj,3,600,20)

    gensc(LMN_list,jj,kk,ll,mm,nn)

    sys.stdout = saved
    fout.close()
