'''General functions for generating atomic structures'''
# from crystal import *
# from defect import *
from gen_pop_box import *
from gen_pop_plate import *
from gen_pop_sphere import *
from generate_dumbbells import *
from get_population import *
from get_restart_population import *
####CHANGE THIS IMPORT TO USE STEM
try:
    from Individual_non_stem import *
except NameError:
    print "NOTE: ASE is not installed. ASE must be installed for Structopt Individual.py to work correctly."
from rot_vec import *
# from surface import *
import crystal
import defect
import surface
