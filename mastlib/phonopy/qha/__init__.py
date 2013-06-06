#!/usr/bin/env python

# Copyright (C) 2012 Atsushi Togo
# All rights reserved.
#
# This file is part of phonopy.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the phonopy project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

try:
    from scipy.optimize import leastsq
except ImportError:
    print "You need to install python-scipy."
    exit(1)

import numpy as np
from phonopy.units import Avogadro, EvTokJmol, EVAngstromToGPa

# Third-order Birch-Murnaghan EOS
def birch_murnaghan(p, v):
    """
    p[0] = E_0
    p[1] = B_0
    p[2] = B'_0
    p[3] = V_0
    """
    return p[0]+9.0/16*p[3]*p[1]*(((p[3]/v)**(2.0/3)-1)**3*p[2]+
                                   ((p[3]/v)**(2.0/3)-1)**2*
                                   (6-4*(p[3]/v)**(2.0/3)))

# Murnaghan EOS
def murnaghan(p, v):
    """
    p[0] = E_0
    p[1] = B_0
    p[2] = B'_0
    p[3] = V_0
    """
    return p[0]+p[1]*v/p[2]*((p[3]/v)**p[2]/(p[2]-1)+1)-p[1]*p[3]/(p[2]-1)

# Vinet EOS
def vinet(p, v):
    """
    p[0] = E_0
    p[1] = B_0
    p[2] = B'_0
    p[3] = V_0
    """
    x = (v / p[3]) ** (1.0 / 3)
    xi = 3.0 / 2 * (p[2] - 1)
    return p[0] + 9 * p[1] * p[3] / (xi**2) * (1 + (xi * (1 - x) - 1) * np.exp(xi * (1 - x)))

def residuals(p, eos, v, e):
    return eos(p, v) - e

class EOSFit:
    def __init__(self, volume, energy, eos=vinet):
        self._energy = np.array(energy)
        self._volume = np.array(volume)
        self._eos = eos

    def run(self, initial_parameter):
        result = leastsq(residuals,
                         initial_parameter[:],
                         args=(self._eos, self._volume, self._energy),
                         full_output=1)
        self._parameters = result[0]

        # self._parameters = [result[0][0],                         # Energy
        #                    result[0][1] * self._pressure_factor,  # Bulk modulus
        #                    result[0][2],                         # Pressure derivative of Bulk modules
        #                    result[0][3]]                        # Volume
        # covar = result[1]
        # print np.sum(self._residuals(result[0], self._volume, self._energy) ** 2)/self._volume.shape[0]
        # print ("%14.6f"*4) % tuple(result[0])
        # print "error              :", ("%14.6f"*4) % tuple(np.sqrt(np.diag(covar)))

    def get_energy(self):
        return self._parameters[0]

    def get_volume(self):
        return self._parameters[3]

    def get_bulk_modulus(self):
        return self._parameters[1]

    def get_parameters(self):
        return self._parameters

def fit_to_eos(volumes, fe, eos=vinet):
    fit = EOSFit(volumes, fe, eos)
    fit.run([fe[len(fe)/2], 1.0, 4.0, volumes[len(volumes)/2]])
    ev = fit.get_volume()
    ee = fit.get_energy()
    eb = fit.get_bulk_modulus()
    ep = fit.get_parameters()
    return ev, ee, eb, ep

def get_bulk_modulus(volumes,
                     electronic_energies,
                     eos=vinet):
    """Returns bulk modulus
    
      (equilibrium volume,
       lowest energy,
       bulk modulus,
       fitting parameters)

    """
    return fit_to_eos(volumes, electronic_energies, eos)    

class QHA:
    def __init__(self,
                 volumes, # Angstrom^3
                 electronic_energies, # eV
                 temperatures,
                 cv,        # J/K/mol
                 entropy,   # J/K/mol
                 fe_phonon, # eV
                 eos=vinet,
                 max_t_index=None,
                 energy_factor=1.0 / EvTokJmol,
                 pressure_factor=EVAngstromToGPa):
        if max_t_index:
            self._max_t_index = max_t_index
        else:
            self._max_t_index = len(temperatures) - 2
        self._volumes = volumes
        self._electronic_energies = electronic_energies
        self._temperatures = temperatures
        self._cv = cv
        self._entropy = entropy
        self._fe_phonon = fe_phonon
        self._eos = eos
        self._energy_factor = energy_factor
        self._pressure_factor = pressure_factor

        self._equiv_volumes = []
        self._equiv_energies = []
        self._equiv_bulk_modulus = []
        self._equiv_parameters = []
        self._free_energies = []

        self._thermal_expansions = None
        self._volume_expansions = None
        self._cp_numerical = None
        self._volume_entropy_parameters = None
        self._volume_cv_parameters = None
        self._volume_entropy = None
        self._volume_cv = None
        self._cp_polyfit = None
        self._dsdv = None
        self._gruneisen_parameters = None
    
    def run(self, verbose=False):
        if verbose:
            print ("#%11s" + "%14s"*4) % ("T", "E_0", "B_0", "B'_0", "V_0")

        for i in range(self._max_t_index + 2):
            t = self._temperatures[i]
            fe = []
            for j, e in enumerate(self._electronic_energies):
                fe.append(e + self._fe_phonon[i][j])
            self._free_energies.append(fe)
    
            ev, ee, eb, ep = fit_to_eos(self._volumes, fe, self._eos)
            self._equiv_volumes.append(ev)
            self._equiv_energies.append(ee)
            self._equiv_bulk_modulus.append(eb * self._pressure_factor)
            self._equiv_parameters.append(ep)

            if verbose:
                print "%14.6f"*5 % (t,
                                    ep[0],
                                    ep[1] * self._pressure_factor,
                                    ep[2],
                                    ep[3])

        self._set_volume_expansion()
        self._set_thermal_expansion() # len = len(t) - 1
        self._set_heat_capacity_P_numerical() # len = len(t) - 2
        self._set_heat_capacity_P_polyfit()
        self._set_gruneisen_parameter() # To be run after thermal expansion.
    
    def plot_helmholtz_volume(self, plt, thin_number=10, display=False):
        if display:
            plt.subplot(1, 3, 1)
        else:
            plt.rcParams['backend'] = 'PDF'
            plt.rcParams['pdf.fonttype'] = 42
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['axes.labelsize'] = 18
            plt.rcParams['figure.subplot.left'] = 0.25
            plt.rcParams['figure.subplot.bottom'] = 0.15
            plt.rcParams['figure.figsize'] = 4, 6

        volume_points = np.linspace(min(self._volumes),
                                    max(self._volumes),
                                    201)
        selected_volumes = []
        selected_energies = []
    
        for i, t in enumerate(self._temperatures[:self._max_t_index]):
            if i % thin_number == 0:
                selected_volumes.append(self._equiv_volumes[i])
                selected_energies.append(self._equiv_energies[i])
                plt.plot(self._volumes,
                         self._free_energies[i],
                         'bo', markersize=4)
                plt.plot(volume_points,
                         self._eos(self._equiv_parameters[i],
                                   volume_points),
                         'b-')

        plt.plot(selected_volumes,
                 selected_energies,
                 'rx-', markersize=8)
        plt.xlabel(r'Volume [$\AA^3$]')
        plt.ylabel('Free energy')
    
        if not display:
            plt.savefig('helmholtz-volume.pdf')
            plt.close()

    def write_helmholtz_volume(self):
        w = open('helmholtz-volume.dat', 'w')
        for i in range(self._max_t_index):
            w.write("# Temperature: %f\n" % self._temperatures[i])
            w.write("# Parameters: %f %f %f %f\n" %
                    tuple(self._equiv_parameters[i]))
            for j, v in enumerate(self._volumes):
                w.write("%20.15f %25.15f\n" % (v, self._free_energies[i][j]))
            w.write("\n\n")
        w.close()
    
    def plot_volume_temperature(self, plt, exp_data=None, display=False):
        if display:
            plt.subplot(1, 3, 2)
        else:    
            plt.rcParams['backend'] = 'PDF'
            plt.rcParams['pdf.fonttype'] = 42
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['axes.labelsize'] = 18
            plt.rcParams['figure.subplot.left'] = 0.15
            plt.rcParams['figure.subplot.bottom'] = 0.15
            plt.rcParams['figure.figsize'] = 8, 6
        
        plt.xlabel('Temperature [K]')
        plt.ylabel(r'Volume [$\AA^3$]')
        plt.plot(self._temperatures[:self._max_t_index],
                 self._equiv_volumes[:self._max_t_index],
                 'b-')
        # exp
        if exp_data:
            plt.plot(exp_data[0], exp_data[1], 'ro')
    
        if not display:
            plt.savefig('volume-temperature.pdf')
            plt.close()
    
    def write_volume_temperature(self):
        w = open('volume-temperature.dat', 'w')
        for i in range(self._max_t_index):
            w.write("%25.15f %25.15f\n" % (self._temperatures[i],
                                           self._equiv_volumes[i]))
        w.close()
    
    def plot_thermal_expansion(self, plt, display=False):
        if display:
            plt.subplot(1, 3, 3)
        else:
            plt.rcParams['backend'] = 'PDF'
            plt.rcParams['pdf.fonttype'] = 42
            plt.rcParams['font.family'] = 'serif'
            plt.rcParams['axes.labelsize'] = 18
            plt.rcParams['figure.subplot.left'] = 0.15
            plt.rcParams['figure.subplot.bottom'] = 0.15
            plt.rcParams['figure.figsize'] = 8, 6
    
        beta = np.array(self._thermal_expansions) * 1e6
        plt.plot(self._temperatures[:self._max_t_index],
                 beta[:self._max_t_index],
                 'b-')
        plt.xlabel('Temperature [K]')
        plt.ylabel('Thermal expansion x 1.0e6 [K^-1]')
    
        if not display:
            plt.savefig('thermal_expansion.pdf')
            plt.close()
    
    def write_thermal_expansion(self):
        w = open('thermal_expansion.dat', 'w')
        for i in range(self._max_t_index):
            w.write("%25.15f %25.15f\n" % (self._temperatures[i],
                                           self._thermal_expansions[i]))
        w.close()
    
    def plot_volume_expansion(self, plt, exp_data=None, symbol='o'):
        plt.rcParams['backend'] = 'PDF'
        plt.rcParams['pdf.fonttype'] = 3
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['axes.labelsize'] = 18
        plt.rcParams['figure.subplot.left'] = 0.15
        plt.rcParams['figure.subplot.bottom'] = 0.15
        plt.rcParams['figure.figsize'] = 8, 6

        plt.plot(self._temperatures[:self._max_t_index],
                 self._volume_expansions[:self._max_t_index],
                 'b-')

        if exp_data:
            plt.plot(exp_data[0],
                     (exp_data[1] / exp_data[1][0]) ** (1.0 / 3) - 1,
                     symbol)

        plt.xlabel('Temperature [K]')
        plt.ylabel(r'Volume expansion $\Delta L/L_0 \, (L=V^{\,1/3})$')
        plt.xlim(self._temperatures[0],
                 self._temperatures[self._max_t_index])
        plt.savefig('volume_expansion.pdf')
        plt.close()
    
    def write_volume_expansion(self):
        w = open('volume_expansion.dat', 'w')
        for i in range(self._max_t_index):
            w.write("%20.15f %25.15f\n" % (self._temperatures[i],
                                           self._volume_expansions[i]))
        w.close()

    def plot_gibbs_temperature(self, plt):
        plt.rcParams['backend'] = 'PDF'
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['axes.labelsize'] = 18
        plt.rcParams['figure.subplot.left'] = 0.15
        plt.rcParams['figure.subplot.bottom'] = 0.15
        plt.rcParams['figure.figsize'] = 8, 6
        plt.xlabel('Temperature [K]')
        plt.ylabel('Gibbs free energy')
        plt.plot(self._temperatures[:self._max_t_index],
                 self._equiv_energies[:self._max_t_index],
                 'b-')
        plt.savefig('gibbs-temperature.pdf')
        plt.close()

    def write_gibbs_temperature(self):
        w = open('gibbs-temperature.dat', 'w')
        for i in range(self._max_t_index):
            w.write("%20.15f %25.15f\n" % (self._temperatures[i],
                                           self._equiv_energies[i]))
        w.close()

    def plot_bulk_modulus_temperature(self, plt):
        plt.rcParams['backend'] = 'PDF'
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['axes.labelsize'] = 18
        plt.rcParams['figure.subplot.left'] = 0.15
        plt.rcParams['figure.subplot.bottom'] = 0.15
        plt.rcParams['figure.figsize'] = 8, 6
        plt.xlabel('Temperature [K]')
        plt.ylabel('Bulk modulus')
        plt.plot(self._temperatures[:self._max_t_index],
                 self._equiv_bulk_modulus[:self._max_t_index],
                 'b-')
        plt.savefig('bulk_modulus-temperature.pdf')
        plt.close()
        
    def write_bulk_modulus_temperature(self):
        w = open('bulk_modulus-temperature.dat', 'w')
        for i in range(self._max_t_index):
            w.write("%20.15f %25.15f\n" % (self._temperatures[i],
                                           self._equiv_bulk_modulus[i]))
        w.close()
    
    def plot_heat_capacity_P_numerical(self, plt, exp_data=None):
        plt.rcParams['backend'] = 'PDF'
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['axes.labelsize'] = 18
        plt.rcParams['figure.subplot.left'] = 0.15
        plt.rcParams['figure.subplot.bottom'] = 0.15
        plt.rcParams['figure.figsize'] = 8, 6
        plt.xlabel('Temperature [K]')
        plt.ylabel(r'$\mathrm{C_P [J/mol\cdot K]}$')
        plt.plot(self._temperatures[:self._max_t_index],
                 self._cp_numerical[:self._max_t_index],
                 'b-')

        # exp
        if exp_data:
            plt.plot(exp_data[0], exp_data[1], 'ro')
                
        plt.savefig('Cp-temperature.pdf')
        plt.close()

    def write_heat_capacity_P_numerical(self):
        w = open('Cp-temperature.dat', 'w')
        for i in range(self._max_t_index):
            w.write("%20.15f %20.15f\n" % (self._temperatures[i],
                                           self._cp_numerical[i]))
        w.close()

    def plot_heat_capacity_P_polyfit(self, plt, exp_data=None):
        plt.rcParams['backend'] = 'PDF'
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['axes.labelsize'] = 18
        plt.rcParams['figure.subplot.left'] = 0.15
        plt.rcParams['figure.subplot.bottom'] = 0.15
        plt.rcParams['figure.figsize'] = 8, 6
        plt.xlabel('Temperature [K]')
        plt.ylabel(r'$\mathrm{C_P [J/mol\cdot K]}$')
        plt.plot(self._temperatures[:self._max_t_index],
                 self._cp_polyfit[:self._max_t_index])
                
        # exp
        if exp_data:
            plt.plot(exp_data[0], exp_data[1], 'ro')

        plt.savefig('Cp-temperature_polyfit.pdf')
        plt.close()

    def write_heat_capacity_P_polyfit(self):
        wve = open('entropy-volume.dat', 'w')
        wvcv = open('Cv-volume.dat', 'w')
        for i in range(1, self._max_t_index):
            t = self._temperatures[i]
            wve.write("# temperature %20.15f\n" % t)
            wve.write("# %20.15f %20.15f %20.15f %20.15f %20.15f\n" %
                      tuple(self._volume_cv_parameters[i - 1]))
            wvcv.write("# temperature %20.15f\n" % t)
            wvcv.write("# %20.15f %20.15f %20.15f %20.15f %20.15f\n" %
                       tuple(self._volume_entropy_parameters[i - 1]))
            for ve, vcv in zip(self._volume_entropy[i - 1],
                               self._volume_cv[i - 1]):
                wve.write("%20.15f %20.15f\n" % tuple(ve))
                wvcv.write("%20.15f %20.15f\n" % tuple(vcv))
            wve.write("\n\n")
            wvcv.write("\n\n")
        wve.close()
        wvcv.close()
    
        w = open('Cp-temperature_polyfit.dat', 'w')
        for i in range(self._max_t_index):
            w.write("%20.15f %20.15f\n" % (self._temperatures[i],
                                           self._cp_polyfit[i]))
        w.close()
    
        w = open('dsdv-temperature.dat', 'w') # GPa
        for i in range(self._max_t_index):
            w.write("%20.15f %20.15f\n" % (self._temperatures[i],
                                           self._dsdv[i] * 1e21 / Avogadro))
        w.close()
    
    def plot_gruneisen_temperature(self, plt):
        plt.rcParams['backend'] = 'PDF'
        plt.rcParams['pdf.fonttype'] = 42
        plt.rcParams['font.family'] = 'serif'
        plt.rcParams['axes.labelsize'] = 18
        plt.rcParams['figure.subplot.left'] = 0.15
        plt.rcParams['figure.subplot.bottom'] = 0.15
        plt.rcParams['figure.figsize'] = 8, 6
        plt.xlabel('Temperature [K]')
        plt.ylabel('Gruneisen parameter')
        plt.plot(self._temperatures[:self._max_t_index],
                 self._gruneisen_parameters[:self._max_t_index],
                 'b-')
        plt.savefig('gruneisen-temperature.pdf')
        plt.close()
        
    def write_gruneisen_temperature(self):
        w = open('gruneisen-temperature.dat', 'w')
        for i in range(self._max_t_index):
            w.write("%20.15f %25.15f\n" % (self._temperatures[i],
                                           self._gruneisen_parameters[i]))
        w.close()
    
    def _set_thermal_expansion(self):
        beta = [0.]
        dt = self._temperatures[1] - self._temperatures[0]
        for i in range(self._max_t_index):
            beta.append((self._equiv_volumes[i + 2] - self._equiv_volumes[i]) /
                        (2 * dt) / self._equiv_volumes[i + 1])

        self._thermal_expansions = beta

    def _set_volume_expansion(self):
        l = np.array(self._equiv_volumes) ** (1.0 / 3)
        for i in range(self._max_t_index):
            t = self._temperatures[i]
            if abs(t - 300) < (self._temperatures[1] - self._temperatures[0]) / 10:
                l_0 = (self._equiv_volumes[i]) ** (1.0 / 3)
                break
    
        self._volume_expansions = l / l_0 - 1

    def _set_heat_capacity_P_numerical(self):
        cp = []
        g = np.array(self._equiv_energies) / self._energy_factor * 1000
        cp.append(0.0)
        cp.append(0.0)
        dt = self._temperatures[1] - self._temperatures[0]
        for i in range(2, self._max_t_index):
            cp.append(-(g[i + 2] - 2 * g[i] + g[i - 2]) /
                       (dt ** 2) / 4 * self._temperatures[i])
        self._cp_numerical = cp
    
    def _set_heat_capacity_P_polyfit(self):
        cp = [0.0]
        dsdv = [0.0]
        wve = open('entropy-volume.dat', 'w')
        wvcv = open('Cv-volume.dat', 'w')
        self._volume_entropy_parameters = []
        self._volume_cv_parameters = []
        self._volume_entropy = []
        self._volume_cv = []
        
        dt = self._temperatures[1] - self._temperatures[0]
        for j in range(1, self._max_t_index):
            t = self._temperatures[j]
            x = self._equiv_volumes[j]
    
            parameters = np.polyfit(self._volumes, self._cv[j], 4)
            cv_p = np.dot(parameters, np.array([x**4, x**3, x**2, x, 1]))
            self._volume_cv_parameters.append(parameters)
    
            parameters = np.polyfit(self._volumes, self._entropy[j], 4)
            dsdv_t = np.dot(parameters[:4], np.array([4*x**3, 3*x**2, 2*x, 1]))
            self._volume_entropy_parameters.append(parameters)
    
            dvdt = (self._equiv_volumes[j + 1] -
                    self._equiv_volumes[j - 1]) / dt / 2
    
            cp.append(cv_p + t * dvdt * dsdv_t)
            dsdv.append(dsdv_t)
    
            self._volume_cv.append(np.array([self._volumes, self._cv[j]]).T)
            self._volume_entropy.append(np.array([self._volumes, self._entropy[j]]).T)

        self._cp_polyfit = cp
        self._dsdv = dsdv

    def _set_gruneisen_parameter(self):
        gamma = [0]
        for i in range(1, self._max_t_index):
            t = self._temperatures[i]
            v = self._equiv_volumes[i]
            kt = self._equiv_bulk_modulus[i]
            beta = self._thermal_expansions[i]
            parameters = np.polyfit(self._volumes, self._cv[i], 4)
            cv = (np.dot(parameters, [v**4, v**3, v**2, v, 1]) /
                  v / 1000 * self._energy_factor * self._pressure_factor)
            gamma.append(beta * kt / cv)
        self._gruneisen_parameters = gamma

