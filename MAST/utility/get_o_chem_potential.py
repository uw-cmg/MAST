#!/usr/bin/env python

######################
# TTM 8/5/15
#     2015-09-24 fix reference mu_O2 to match MP
#     2015-10-23 add in VASP energy terms
#######################

#####################
# References
# http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1#Thermo-Gas
# Chase, M.W., Jr., NIST-JANAF Themochemical Tables, Fourth Edition, J. Phys. Chem. Ref. Data, Monograph 9, 1998, 1-1951
# Cox, J.D.; Wagman, D.D.; Medvedev, V.A., CODATA Key Values for Thermodynamics, Hemisphere Publishing Corp., New York, 1984, 1. 
# Ping Ong, S., Wang, L., Kang, B. & Ceder, G. Li-Fe-P-O2 Phase Diagram from First Principles Calculations. Chem Mater 20, 1798-1807, doi:10.1021/cm702327g (2008).
# Ong, S. P., Jain, A., Hautier, G., Kang, B. & Ceder, G. Thermal stabilities of delithiated olivine MPO4 (M=Fe, Mn) cathodes investigated using first principles calculations. Electrochemistry Communications 12, 427-430, doi:10.1016/j.elecom.2010.01.010 (2010).

######################
# S0 gas, 1 bar 205.15(2) J/molK
#####################
# 0.1 MPa = 1 bar
#####################

import sys

##################
# Arguments: 
# temperature in Kelvin
# oxygen partial pressure in atm
################

if len(sys.argv) < 3:
    print "Run with the first parameter as temperature in K"
    print "and the second parameter as oxygen partial pressure in atm."
    print "Example: python get_o_chem_potential.py 1073 0.21"
    sys.exit()

temp_K = float(sys.argv[1])
pressure_atm = float(sys.argv[2])

#T=100-700K
dict1=dict()
dict1["A"]= 31.32234
dict1["B"]= -20.23531
dict1["C"]= 57.86644
dict1["D"]= -36.50624
dict1["E"]= -0.007374
dict1["F"]= -8.903471
dict1["G"]= 246.7945 
dict1["H"]= 0

#T=700-2000K
dict2=dict()
dict2["A"]=    30.03235
dict2["B"]=   8.772972 
dict2["C"]=    -3.988133
dict2["D"]=   0.788313
dict2["E"]=   -0.741599
dict2["F"]=   -11.32468
dict2["G"]=    236.1663
dict2["H"]=   0  

#T=2000-6000K
dict3=dict()
dict3["A"]=    20.91111
dict3["B"]=    10.72071
dict3["C"]=   -2.020498
dict3["D"]=    0.146449
dict3["E"]=   9.245722
dict3["F"]=   5.337651
dict3["G"]=   237.6185
dict3["H"]=  0
             

if temp_K < 700:
    cdict=dict(dict1)
elif temp_K < 2000:
    cdict=dict(dict2)
else:
    cdict=dict(dict3)

t = temp_K/1000.0
t2 = t*t
t3 = t*t*t
t4 = t*t*t*t

#heat capacity (J/mol*K)
#Cp0=A + B*t + C*t2 + D*t3 + E/t2
heat_capacity = cdict["A"] + cdict["B"]*t + cdict["C"]*t2 + cdict["D"]*t3 + cdict["E"]/t2

#standard enthalpy (kJ/mol)
#H0-H0_298.15=A*t + B*t2/2 + C*t3/3 + D*t4/4 - E/t + F - H
enthalpy=cdict["A"]*t + cdict["B"]*t2/2.0 + cdict["C"]*t3/3.0 + cdict["D"]*t4/4.0 - cdict["E"]/t + cdict["F"] - cdict["H"]

#standard entropy (J/mol*K)
#S0=A*ln(t) + B*t + C*t2/2 + D*t3/3 - E/(2*t2) + G
import numpy as np
entropy=cdict["A"] * np.log(t) + cdict["B"]*t + cdict["C"]*t2/2.0 + cdict["D"]*t3/3.0 - cdict["E"]/(2*t2) + cdict["G"]


kJpermol_in_eVperatom=96.487
Cp_eVperK = heat_capacity/1000.0/kJpermol_in_eVperatom
H_eV = enthalpy/kJpermol_in_eVperatom
S_eVperK = entropy/1000.0/kJpermol_in_eVperatom

print "Cp (J/mol*K)", heat_capacity, " (eV/atomK): ", Cp_eVperK
print "H0-H0_298.15 (kJ/mol)", enthalpy, " (eV/atom): ", H_eV
print "S0 (J/mol*K)", entropy, " (eV/atomK): ", S_eVperK

pascals_in_atm = 101325
pressure_Pa = pressure_atm*101325.0
pressure_ref_Pa = 100000.0 #1 bar, equivalent to 0.1 MPa, for NIST data

boltzmann = 0.00008617

#muO2(T,pO2)=muO2(T,p0) + kTln(pO2/p0)
#=approx= h_O2(T,p0) -T(s_O2(T,p0)-kln(pO2/p0))



#E_O2_VASP = -9.0903605 # eV, PAW-PW91 O_s, Morgan group open systems document
E_O2_VASP = -9.4027177 # eV, PAW-PBE O_s, Morgan group open systems document
#E_O2_VASP = -9.7819830 # eV, PAW-PW91 O, Morgan group open systems document
#E_O2_VASP = -9.8644451 # eV, PAW-PBE O, Morgan group open systems document

E_O2_VASP_shift = 1.36 # eV, Wang et al, Physical Review B 73, 195107 (2006)
#In the Morgan group open systems document, the shift is applied as +1.36 eV.
#However, from Wang (2006) it seems like the shift should be -1.36 eV.
#Shifting up is correct to correct for GGA tendency towards reducing conditions (so the correction should make solids more stable and O2 gas less stable)

E_O2_VASP_corrected = E_O2_VASP + E_O2_VASP_shift


first_term = E_O2_VASP_corrected # T=0 energy of O2
second_term = H_eV - temp_K * S_eVperK #relative free energy of O2 gas relative to gas enthalpy at ref_pressure and ref_temperature
#The H_eV term already has enthalpy relative to ref_pressure and ref_temperature subtracted.

third_term = temp_K * boltzmann * np.log(pressure_Pa/pressure_ref_Pa) #pressure effect


mu_O2 = first_term + second_term + third_term
mu_O = 0.5 * mu_O2

#mpadjust = 8.716 #adjustment for Materials Project so that mu_O2 reference is 0 at 298K and 0.21 atm pO2 for Anubhav's 2010 paper
#mu_O2_adjusted = first_term + second_term + mpadjust + third_term
#mu_O_adjusted = 0.5 * mu_O2_adjusted

print "mu_O2 (eV): ", mu_O2
print "mu_O (eV): ", mu_O

