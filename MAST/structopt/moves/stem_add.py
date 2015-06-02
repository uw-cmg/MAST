import random
from MAST.structopt.tools.StemCalc import ConvStem
from ase import Atoms, Atom
import math
from math import fabs, sqrt, acos, pi, atan
import numpy 

def On_boundary(map_atom_pixel,x_pixel,y_pixel):
    
    for x in [-2,-1,0,1,2]:
       for y in [-2,-1,0,1,2]:
          if map_atom_pixel[x_pixel+x][y_pixel+y] > 5: 
            return False
    return True        

def stem_add(indiv, Optimizer):
    atms = indiv[0].copy()
    simfun = Optimizer.stemcalc.get_image(Optimizer.stemcalc.psf,atms,Optimizer.stemcalc.parameters['Slice size'],Optimizer.stemcalc.parameters['Pixels'])
    expfun = Optimizer.stemcalc.expfun

    xmax = Optimizer.stemcalc.parameters['Slice size']
    com = atms.get_center_of_mass()
    #com = [ 44.40963074 , 44.65497562 , 44.90406073]
    #com = numpy.array(com)
    #com += [-0.149836425, 0.29967285, 0]  #pixelshift
    atm_num = atms.get_number_of_atoms()
    cop = xmax/2.0 
    trans = [cop-i for i in com]
    atms.translate(trans)
    R = atms.arrays['positions']
    ax=[]
    ay=[]
    az=[] 
    for o,t,h in R:
        ax.append(o)
        ay.append(t)
        az.append(h)
    ax = numpy.array(ax)
    ay = numpy.array(ay)
    az = numpy.array(az)

    nx = Optimizer.stemcalc.parameters['Pixels']
    ny = nx
    dx = xmax/nx
    dy = dx
    ix = numpy.array([math.floor(axi/dx) for axi in ax])
    iax = numpy.array([int(math.fmod(iaxi,nx)) for iaxi in ix])
    ibx = numpy.array([int(math.fmod(iaxi+1,nx)) for iaxi in ix])
    iy = numpy.array([math.floor(ayi/dy) for ayi in ay])
    iay = numpy.array([int(math.fmod(iayi,ny)) for iayi in iy])
    iby = numpy.array([int(math.fmod(iayi+1,ny)) for iayi in iy])

 
# construct a map from atom to pixel grid   
    map_atom_pixel = []

    for i in range(nx) :
       map_atom_pixel.append([])
       for j in range(ny) :
          map_atom_pixel[i].append([])

    for iatom in range(len(ix)):
       for x in [iax[iatom],ibx[iatom]] :
          for y in [iay[iatom],iby[iatom]] :
              map_atom_pixel[x][y].append(iatom)
 
#add atoms to negative delta STEM intensity 
    delta = simfun-expfun
    for i in range(len(expfun)):
        for j in range(len(expfun[0])):
           if expfun[i][j] <= -1000 : #invalid exp points
             delta[i][j] = 0 
    add_tol = delta.min()*0.3
    maxmove = 1

# construct the pixel list where simfun differs from expfun
    add_list_x = []
    add_list_y = []

    for x in range(nx) :
       for y in range(ny) :
          if delta[x][y] < add_tol:
             add_list_x.append(x)
             add_list_y.append(y)
    add_list_x = numpy.array(add_list_x)
    add_list_y = numpy.array(add_list_y)

#adding atom to delta negative col 
    addcollist = []
    for imove in range(maxmove): 
       addcol = random.randint(0,len(add_list_x)-1)
       addcollist.append(addcol)
    addcollist = list(set(addcollist))
    for addcol in addcollist:     
       x_pixel = add_list_x[addcol]
       y_pixel = add_list_y[addcol]
       #print 'add',x_pixel,y_pixel,delta[x_pixel][y_pixel],map_atom_pixel[x_pixel][y_pixel]
       #print 'add',numpy.array([R[j][2] for j in map_atom_pixel[x_pixel][y_pixel]])
       if len(map_atom_pixel[x_pixel][y_pixel]) > 0:
         ind = map_atom_pixel[x_pixel][y_pixel][0]
         atomind_max = ind
         atomind_min = ind
         Zmax = R[ind][2]-cop
         Zmin = R[ind][2]-cop
         for j in range(1,len(map_atom_pixel[x_pixel][y_pixel])):
           ind = map_atom_pixel[x_pixel][y_pixel][j]
           if R[ind][2]-cop > Zmax :
             atomind_max = ind
             Zmax = R[ind][2]-cop
           elif R[ind][2]-cop < Zmin : 
             atomind_min = ind
             Zmin = R[ind][2]-cop  
         atomind_add = atomind_max #random.choice([atomind_max,atomind_min]) 
         #ignore this addcol if one atom has been added to this col already.
         if atomind_add >= atm_num :
            continue
         if atomind_add == atomind_max :
            sign = 1
         elif atomind_add == atomind_min :
            sign = -1
         pos = [R[atomind_add][0],R[atomind_add][1],R[atomind_add][2]+sign*2.6]
         atomsym=random.choice(Optimizer.atomlist)[0]
         at=Atom(atomsym,pos)     
         atms.append(at)
         R = atms.arrays['positions']
         #print 'add',atomind_add,R[atomind_add][2],pos
       # remove this atom index from removecol of map_atom_pixel  
         for x in [iax[atomind_add],ibx[atomind_add]] :
            for y in [iay[atomind_add],iby[atomind_add]] :
               map_atom_pixel[x][y].append(atms.get_number_of_atoms()-1)
       elif On_boundary(map_atom_pixel,x_pixel,y_pixel) :
         pos = [x_pixel*dx,y_pixel*dy,cop]
         atomsym=random.choice(Optimizer.atomlist)[0]
         at=Atom(atomsym,pos)     
         atms.append(at)
         R = atms.arrays['positions']
         #print 'add',pos
         map_atom_pixel[x_pixel][y_pixel].append(atms.get_number_of_atoms()-1) 
        
 
    indiv[0] = atms
    muttype='STEM'+repr(maxmove)
    Optimizer.output.write('STEM add performed on individual\n')
    #Optimizer.output.write('Index = '+repr(add)+'\n')
    #Optimizer.output.write(repr(indiv[0])+'\n')
    #Optimizer.output.write(repr(indiv[0].get_positions())+'\n')
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype

    return indiv
