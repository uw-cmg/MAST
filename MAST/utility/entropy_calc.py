import sys
import math

#print 'Number of arguments:', len(sys.argv), 'arguments.'

def get_O2_entropy(temp,pressure):
    #Modified from Brian Puchala's MATLAB code:
    
    #Constants
    #All parameters from NIST webbook: http://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Units=SI&Type=JANAFG&Table=on#JANAFG
    
    #Boltzmann constant [eV/K]
    k_b = 8.617343e-5
    
    #1eV = 96486 [J/mole]
    eV = 96486
    
    #1 bar = 0.98692 [atm]
    bar = 0.98692
    
    #VASP O2 energy (PAW-PBE O) [eV]
    E_O2 = -9.8644451
    
    #Correction of O2 energy [eV], from Wang et al, PHYSICAL REVIEW B 73, 195107 (2006)
    E_O2_corr = 1.36
    
    temp_mod = temp/1000
    
    if temp < 100:
        print 'Low Temperature!'
        return
    if temp >= 100 and temp < 700:
        A = 31.32234
        B = -20.23531
        C = 57.86644
        D = -36.50624
        E = -0.007374
        F = -8.903471
        G = 246.7945
    if temp >= 700 and temp < 2000:
        A = 30.03235
        B = 8.772972
        C = -3.988133
        D = 0.788313
        E = -0.741599
        F = -11.32468
        G = 236.1663
    if temp >= 2000 and temp <= 6000:
        A = 20.91111
        B = 10.72071
        C = -2.020498
        D = 0.146449
        E = 9.245722
        F = 5.337651
        G = 237.6185
    if temp > 6000:
        print 'High Temperature!'
        return
    #Entropy of oxygen gas (J/mol*K, 298~6000K): Shomate Equation
    S_O2 = ( A*math.log(temp_mod) + B*temp_mod + C*math.pow(temp_mod,2)/2 + D*math.pow(temp_mod,3)/3 - (E/(2*math.pow(temp_mod,2))) + G )/eV
    S_O2_std = 205.0694/eV # S_O2 at temp = 298 K
# if have a temp between 296 K and 300 K, can see that S_O2_std is S_O2 at 298 K:
#    if temp >= 296 and temp <= 300:
#        print temp
#        print S_O2
#        print S_O2_std
    
    #(reference to H_O2_298.15)
    H_O2_NIST = ( A*temp_mod + B*math.pow(temp_mod,2)/2 + C*math.pow(temp_mod,3)/3 + D*math.pow(temp_mod,4)/4 - E/temp_mod + F )*1000/eV

    # Gibbs energy of oxygen gas (G_O2 references to G_O2_298.15)
    G_O2 = ( H_O2_NIST - temp*S_O2 + 298.15*S_O2_std ) + k_b*temp*math.log(pressure/bar)
    
    #Enthalpy of oxygen gas using G=H-TS to obtain H from G.
    H_O2 = ( E_O2 + E_O2_corr ) + G_O2 + ( -298.15*S_O2_std ) + temp*S_O2
    
    #A plain way to get H_O2 is using NIST H_O2 equation and VASP E_O2 plus s correction, which will give the same result. H_O2=(E_O2+E_O2_corr)+ H_O2_NIST
    #print S_O2
    
#    print 'Temp = ', temp, ' [K], Pres = ', pressure, ' [bar], Entropy = ', S_O2, '[eV].'
    print 'Temp = ', temp, ' [K], Entropy = ', S_O2, '[eV].'  
    #return S_O2

#get_O2_entropy(1000,1)
steps = 300
min_temp = 200.
max_temp = 1400.
step_size_temp = ( max_temp - min_temp )/steps
min_press = 1e-9
max_press = 1e-5
step_size_press = ( max_press - min_press )/steps
temp_list = []
press_list = []

for i in xrange(0,steps+1):
    temp_list.append(min_temp + step_size_temp*i)
    press_list.append(min_press + step_size_press*i)

for temp in temp_list:
    pressure = 1e-5
#    for pressure in press_list:
    get_O2_entropy(temp,pressure)
