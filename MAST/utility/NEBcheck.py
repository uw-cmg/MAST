##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Zhewen Song
# Last updated: 2014-04-25
##############################################################
import pymatgen as mg
from numpy import *
import create_paths
import os, getopt, sys, re

argv = sys.argv[1:]
opts, args = getopt.getopt(argv,"i","n")
inp = os.path.join(os.getcwd(),args[0])
Nth = int(args[2])
info=create_paths.main(inp,Nth)
[vec,xyz,path]=info

def add_structure():
    fp=open('temp.inp','w+')
    rp=open(inp,'r+').readlines()
    i=0
    while not 'begin lattice' in rp[i]: i=i+1
    while not 'end' in rp[i+1]: del rp[i+1]
    for j in range(3):
        i=i+1; rp.insert(i,str(vec[j][0])+' '+str(vec[j][1])+' '+str(vec[j][2])+' \n')   
    while not 'begin coordinates' in rp[i]: i=i+1
    while not 'end' in rp[i+1]: del rp[i+1]
    for ele in xyz.keys():
        for j in range(len(xyz[ele])):
            i=i+1; rp.insert(i,ele+' '+str(xyz[ele][j][0])+' '+str(xyz[ele][j][1])+' '+str(xyz[ele][j][2])+'\n')
    while not '$site' in rp[i]: i=i+1
    del rp[i]
    while not '$end' in rp[i]: del rp[i]
    fp.writelines(rp)
    
def add_defects():
    ends=[]
    tmp=[]
    dash=re.compile(r'-')
    digit=re.compile(r'\d')
    point=re.compile(r'\.')
    for i in range(len(path['T'])):
        for j in [2,3]:
            if not path['T'][i][j] in ends:
                ends.append(path['T'][i][j])
                name=dash.split(path['T'][i][0])[j-2]
                tmp.append(name)
    types=tmp[:]
    types[0]=tmp[0]+chr(ord('a'))
    for i in range(1,len(tmp)):
        flag=0
        for j in range(i-1,-1,-1):
            if tmp[j]==types[i]:
                types[i]=tmp[i]+chr(ord(digit.split(types[j])[1])+1)
                flag=1
                break
        if flag==0: 
            types[i]=tmp[i]+chr(ord('a')) # New type
    name=point.split(inp)[0]+'_new.inp'    
    fp=open(name,'w+')
    rp=open('temp.inp','r+').readlines()
    rp.append('$defects\nthreshold 1e-4\ncoord_type fractional\n')
    for i in range(len(ends)):   
        rp.append('begin '+types[i]+'\n')
        rp.append('interstitial '+str(ends[i][0])+' '+str(ends[i][1])+' '+str(ends[i][2])+' '+'X'+str(len(xyz.keys())+1)+'\n')
        rp.append('end\n')
    rp.append('$end\n\n')
    rp.append('$neb\n')
    for i in range(len(path['T'])):
        for j in range(len(ends)):
            if ends[j]==path['T'][i][2]:
                start=j
            if ends[j]==path['T'][i][3]:   
                finish=j
        rp.append('begin '+types[start]+'-'+types[finish]+' #'+str(path['T'][i][1])+'NN\n')
        rp.append('X'+str(len(xyz.keys())+1)+', '+str(ends[start][0])+' '+str(ends[start][1])+' '+str(ends[start][2])+', '+str(ends[finish][0])+' '+str(ends[finish][1])+' '+str(ends[finish][2])+'\n')
        rp.append('images 1\n')
        mid=[0,0,0]
        for k in range(3): 
            mid[k]=(ends[start][k]+ends[finish][k])*0.5
            if abs(ends[start][k]-ends[finish][k])>0.5:
                mid[k]=mid[k]-0.5
            
        rp.append('phonon int '+str(mid[0])+' '+str(mid[1])+' '+str(mid[2])+' 0.0\n')
        rp.append('end\n')
    rp.append('$end\n')
    os.system('rm temp.inp')  
    fp.writelines(rp)
    return name

add_structure()
name=add_defects()
os.system('mast -i '+name)
