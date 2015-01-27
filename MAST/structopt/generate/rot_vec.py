import numpy as np
import math
from math import cos as cos
from math import sin as sin
from math import pi as pi
import random

def rot_vec(pos_list, al, be, ga):
   print 'HKK:: you ar in rot_vec' 
#   '''
#   vec=np.array([0, 1, 0])
#   al=pi/2
#   be=0
#   ga=0
#   '''
   vec=np.asarray(pos_list)
   al=-al
   be=-be
   ga=-ga
   rot_mat=np.array([[cos(be)*cos(ga), cos(al)*sin(ga)+sin(al)*sin(be)*cos(ga),sin(al)*sin(ga)-cos(al)*sin(be)*cos(ga)], [-cos(be)*sin(ga), cos(al)*cos(ga)-sin(al)*sin(be)*sin(ga),sin(al)*cos(ga)+cos(al)*sin(be)*sin(ga)],[sin(be), -sin(al)*cos(be),cos(al)*cos(be)]])
   rotvec=np.dot(rot_mat,vec)
   rotvec=rotvec.tolist()
   print 'final', rotvec
   return rotvec
