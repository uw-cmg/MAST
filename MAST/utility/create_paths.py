##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Zhewen Song
# Last updated: 2014-04-25
##############################################################
from numpy import *
import pymatgen as mg

def getinfo(line):
    line=line.strip('\n')
    data=line.split(' ')
    while 1:
        try: data.remove('')
        except: break
    return data

def read_ele(inp):
    fp=open(inp,'r+').readlines()
    i=0
    ele={}
    while i<len(fp):
        if 'begin elementmap' in fp[i]:
            i=i+1
            while not 'end' in fp[i]:
                line=getinfo(fp[i])
                if len(line)==2:
                    key=line[0]
                    ele[key]=line[1]
                i=i+1
        i=i+1
    return ele
            
            
def read_site(inp):
    fp=open(inp,'r+').readlines()
    i=0
    site={}
    while i<len(fp):
        if '$site' in fp[i]:
            i=i+1
            while not '$end' in fp[i]:
                line=getinfo(fp[i])
                if len(line)==1:
                    key=line[0]
                    site[key]=[]
                elif len(line)==3:
                    site[key].append([float(line[0]),float(line[1]),float(line[2])]) 
                i=i+1
        i=i+1
    return site
    
def read_vec(inp):
    fp=open(inp,'r+').readlines()
    i=0           
    lattice_vectors=zeros((3,3)) 
    while i<len(fp):
        if 'begin lattice' in fp[i]:
            i=i+1; j=0
            while not 'end' in fp[i]:
                line=getinfo(fp[i])
                if len(line)==3:
                    for k in range(3):
                        lattice_vectors[j][k]=float(line[k]) 
                    j=j+1
                i=i+1     
        i=i+1
    return  lattice_vectors
    
def read_coordinates(inp):
    fp=open(inp,'r+').readlines()
    i=0           
    coordinates={}   
    while i<len(fp):
        if 'begin coordinates' in fp[i]:
            i=i+1
            while not 'end' in fp[i]:
                line=getinfo(fp[i])
                if len(line)==4:
                    try: coordinates[line[0]]
                    except: coordinates[line[0]]=[]
                    coordinates[line[0]].append([float(line[1]),float(line[2]),float(line[3])])
                i=i+1
        i=i+1
    return coordinates
        
def MakeSupercell(xyz,size):
    supercell={}
    for i in range(size[0]):
        for j in range(size[1]):
            for k in range(size[2]):
                for ele in xyz.keys():
                    for atom in range(len(xyz[ele])):
                        try: supercell[ele]
                        except: supercell[ele]=[]
                        new=[1.*(xyz[ele][atom][0]+i)/size[0],1.*(xyz[ele][atom][1]+j)/size[1],1.*(xyz[ele][atom][2]+k)/size[2]]
                        supercell[ele].append(new)
    return supercell
       
def distance(A,B,vec):
    d=array([0.,0.,0.])
    for i in range(3):
        d[i]=A[i]-B[i]
        if d[i]>0.5: d[i]=d[i]-1
        elif d[i]<-0.5: d[i]=d[i]+1
    d=dot(d,vec)
    return double(linalg.norm(d))

def neighbors(site,Nth,vec):  # Get nth neighbor of possible pairs of same/different site types
    pair={}
    for i in site.keys():
        for j in site.keys():
            val=[]
            for ii in range(len(site[i])):
                for jj in range(len(site[j])):
                    d=distance(site[i][ii],site[j][jj],vec)
                    if d>1e-5:
                        val.append((d,site[i][ii],site[j][jj]))
            dist=array(val,dtype=[('length',float),('coord1',list),('coord2',list)])  
            if len(dist)==0: return 'NaN'
            dist=sort(dist,order='length') 
            count=1
            k=0
            for k in range(len(dist)-1):
                if dist[k+1][0]-dist[k][0]>1e-3:
                    if count==Nth: 
                        if not j+'-'+i in pair.keys():
                            pair[i+'-'+j]=[dist[k][1],dist[k][2]]
                    count=count+1
            if count==Nth: 
                if not j+'-'+i in pair.keys():
                    pair[i+'-'+j]=[dist[k][1],dist[k][2]]            
            if count<Nth: return 'NaN' # Input cell too small, Nth nearest neighbor not found
    return pair            

def isline(A,M,B,vec):
    a=distance(A,M,vec)
    b=distance(B,M,vec)
    c=distance(A,B,vec)
    p=(a+b+c)/2.0
    area=(p*(p-a)*(p-b)*(p-c))**0.5
    if c>a and c>b:
	if 2*area/c<1e-3: flag=1
        elif 2*area/c<1.0: flag=-1
	else: flag=0
    else: flag=0
    return flag
    
def iscross(xyz,vec,site,site1,site2):  
    for ele in xyz.keys():
        for i in range(len(xyz[ele])): # check if the path crosses over a host atom
            coords=xyz[ele][i]
            if not distance(coords,site1,vec)==0 and not distance(coords,site2,vec)==0:
                flag=isline(site1,coords,site2,vec)
		if not flag==0: return 'a host lattice!'
    for i in site.keys(): # check if the path crosses over another site
        for j in range(len(site[i])):
            if not distance(site[i][j],site1,vec)==0 and not distance(site[i][j],site2,vec)==0: 
                flag=isline(site1,site[i][j],site2,vec)
		if flag==1: return 'site '+i+'!'    
    return 1
            
def get_path(site,xyz,vec,N): 
    path={'T':[],'F':[]}
    for nb in range(1,N+1):
        if neighbors(site,nb,vec)=='NaN':
            return 'NaN'
        else: pair=neighbors(site,nb,vec)
        for i in pair.keys():
            cross=iscross(xyz,vec,site,pair[i][0],pair[i][1])
            if cross==1:
                path['T'].append([i,nb,pair[i][0],pair[i][1]])
            else: 
                path['F'].append([i,nb,pair[i][0],pair[i][1],cross])
    return path
    

def main(inp,N):
    size=array([2,2,2])
    vec=read_vec(inp)
    xyz=read_coordinates(inp)
    site=read_site(inp)
    xyz=MakeSupercell(xyz,size)
    site=MakeSupercell(site,size)
    vec=vec*size
    path=get_path(site,xyz,vec,N) 
    i=1
    while path=='NaN':
        i=i+1
        size=size/(i-1)*i
        xyz=MakeSupercell(xyz,size)
        site=MakeSupercell(site,size)
        vec=vec*size
        path=get_path(site,xyz,vec,N)         
    if False in (size==array([1,1,1])):
        print 'The input cell size is too small to find all '+str(N)+' NN, so the size is expanded to '+str(size[0])+'x'+str(size[1])+'x'+str(size[2])+' supercell!'
    
    print 'There are '+str(len(path['T'])+len(path['F']))+' paths in total.'
    if not len(path['F'])==0:
        print 'But only '+str(len(path['T']))+' paths are possible.'
        for i in range(len(path['F'])):
            print 'Path '+path['F'][i][0]+'('+str(path['F'][i][1])+' NN) is removed because it is too close to '+path['F'][i][4]
     #   print 'These removed paths are labeled as POSCAR_F*.'
    return [vec,xyz,path]
