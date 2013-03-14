#-------------------------------------------------------------------------------
# Name:        Structure
# Purpose:
#
# Author:      Kumaresh Visakan Murugan
#
# Created:     11/01/2013
#-------------------------------------------------------------------------------

import numpy as np
from atom import Atom
import PhysicalConstants as PC
from gaussj import gaussJordan
import useful_funcs as uf

#number of elements
NUM_ELEMENTS    = 110
#elements ordered by atomic number
#F-valence electron atoms and atomic numbers > 86 not found here
AT_NAMES        = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc",\
                   "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb",\
                   "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Hf", "Ta", "W", "Re",\
                   "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn" ]

class Structure:
    #constructor
    def __init__(self):
        self.num_types            = 0
        self.num_atoms            = 0
        self.title                = "NO_TITLE_GIVEN"
        self.scale                = 1
        self.lat                  = np.matrix([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
        self.num_each_type        = np.array([])
        self.types                = np.array([])
        self.coord_type           = ['D', 'D']
        # 0 => direct, 1 => cart
        self.coord_flag           = 0
        self.cpos                 = np.matrix([[]])
        self.dpos                 = np.matrix([[]])
        self.names                = np.array([])
        self.names_were_given     = np.array([])
        self.at_num_vec           = np.array([0 for i in xrange(NUM_ELEMENTS)])
        self.at_name_vec          = np.array(["" for i in xrange(NUM_ELEMENTS)])
        self.initialize_atomic_numbers()
        self.xray_scatt_vec       = np.array([i * 1.0 for i in xrange(NUM_ELEMENTS)])
        self.initialize_xray_scatt_vec()
        self.at_mass_vec          = np.array([2.0 * (i + 1) * PC.amu2kg for i in xrange(NUM_ELEMENTS)])
        # 0 => no selective dynamics, !0 => selective dynamics
        self.isd                  = 0
        #selective dynamic settings
        self.sd                   = np.array([])
        # 0 => scale >= 0, !0 => scale < 0
        self.neg_scale            = 0
        self.origin               = np.array([0.0, 0.0, 0.0])

    def initialize_atomic_numbers(self):
        at_names_len = len(AT_NAMES)
        for i in xrange(NUM_ELEMENTS):
            self.at_num_vec[i]  = i + 1
            if i < at_names_len:
                self.at_name_vec[i] = AT_NAMES[i]


    def initialize_xray_scatt_vec(self):
        #All indices are the atomic number shifted back by one.
        #All data collected from the NIST online tables:
        #http://physics.nist.gov/PhysRefData/FFast/html/form.html
        #All data are ideally for f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
        #These are for E=7.9026keV (Cu-alpha is wavelength=1.5418A, E=8.0416keV).
        self.xray_scatt_vec[2]=3.00145E+00; # Li
        self.xray_scatt_vec[14]=1.53133E+01; # P
        self.xray_scatt_vec[24]=2.43589E+01; # Mn
        self.xray_scatt_vec[25]=2.46830E+01; # Fe

        # All data collected from the online tables:
        # http:#www-cxro.lbl.gov/optical_constants/pert_form.html
        # All data are f1 values for Cu-alpha (wavelength=1.5418A, E=8.0416keV).
        self.xray_scatt_vec[0]=1.000; # H
        self.xray_scatt_vec[1]=2.000; # He
        self.xray_scatt_vec[2]=3.001; # Li
        self.xray_scatt_vec[6]=6.019; # C
        self.xray_scatt_vec[7]=8.052; # O
        self.xray_scatt_vec[13]=14.43; # P
        self.xray_scatt_vec[14]=15.30; # P
        self.xray_scatt_vec[20]=21.34; # Sc
        self.xray_scatt_vec[21]=22.24; # Ti
        self.xray_scatt_vec[23]=23.84; # Cr
        self.xray_scatt_vec[24]=24.46; # Mn
        self.xray_scatt_vec[25]=24.85; # Fe
        self.xray_scatt_vec[26]=24.59; # Co
        self.xray_scatt_vec[27]=25.02; # Ni
        self.xray_scatt_vec[28]=27.03; # Cu
        self.xray_scatt_vec[29]=28.44; # Zn
        self.xray_scatt_vec[46]=47.18; # Ag
        self.xray_scatt_vec[78]=74.99; # Au
        self.xray_scatt_vec[89]=86.64; # Th

    def setTitle(self, in_title):
        self.title = title

    def setScale(self, in_scale):
        self.scale = scale

    def setLat(self, lat):
        self.lat = lat

    def setNumEachType(self, num_each_type):
        self.num_each_type = num_each_type

    def setTypes(self, types):
        self.types = types

    def setCoordType(self, coord_type):
        self.coord_type = coord_type
        if self.coord_type in ["D", "d"]:
            self.coord_flag = 0
        elif self.coord_type in ["C", "c"]:
            self.coord_flag = 1

    def setCoordFlag(self, coord_flag):
        self.coord_flag = coord_flag
        if self.coord_flag == 0:
            self.coord_type = "D"
        elif self.coord_flag == 1:
            self.coord_type = "C"

    def setAllAtomPos(self, pos, coord_flag):
        if coord_flag not in [0, 1]:
            return
        if coord_flag == 0:
            self.dpos      = pos
            self.num_atoms = len(self.dpos)
            self.cpos      = np.matrix([[0.0, 0.0, 0.0] for i in xrange(self.num_atoms)])
            self.calcAllCartCoords()
        elif coord_flag == 1:
            self.c_pos     = pos
            self.num_atoms = len(self.cpos)
            self.dpos      = np.matrix([[0.0, 0.0, 0.0] for i in xrange(self.num_atoms)])
            self.calcAllDirectCoords()

    def rescale(self, in_scale):
        if in_scale == 0:
            return
        for i in xrange(3):
            for j in xrange(3):
                self.lat[i][j] *= self.scale / in_scale

        for i in xrange(self.num_atoms):
            for j in xrange(3):
                self.cpos[i][j] *= self.scale / in_scale

        self.scale = in_scale


    def _str_abs3D(self, x):
        return np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) * 1.0

    def _str_scl3D(self, x, y):
        return ( x[0]*y[0] + x[1]*y[1] + x[2]*y[2] ) * 1.0

    def WignerSeitz():
        rat = [0.0 for i in xrange(3)]
        xoo = [0.0 for i in xrange(3)]
        yoo = [0.0 for i in xrange(3)]
        zoo = [0.0 for i in xrange(3)]
        xyo = [0.0 for i in xrange(3)]
        xzo = [0.0 for i in xrange(3)]
        yzo = [0.0 for i in xrange(3)]
        xyz = [0.0 for i in xrange(3)]

        xoo[0] = self.lat[0][0]
        xoo[1] = self.lat[0][1]
        xoo[2] = self.lat[0][2]
        axoo   = self._str_abs3D(xoo)

        yoo[0] = self.lat[1][0]
        yoo[1] = self.lat[1][1]
        yoo[2] = self.lat[1][2]
        ayoo   = self._str_abs3D(yoo)

        zoo[0] = self.lat[2][0]
        zoo[1] = self.lat[2][1]
        zoo[2] = self.lat[2][2]
        azoo   = self._str_abs3D(zoo)

        xyo[0] = self.lat[0][0] + self.lat[1][0]
        xyo[1] = self.lat[0][1] + self.lat[1][1]
        xyo[2] = self.lat[0][2] + self.lat[1][2]
        axyo   = self._str_abs3D(xyo)

        xzo[0] = self.lat[0][0] + self.lat[2][0]
        xzo[1] = self.lat[0][1] + self.lat[2][1]
        xzo[2] = self.lat[0][2] + self.lat[2][2]
        axyo   = self._str_abs3D(xzo)

        yzo[0] = self.lat[1][0] + self.lat[2][0]
        yzo[1] = self.lat[1][1] + self.lat[2][1]
        yzo[2] = self.lat[1][2] + self.lat[2][2]
        ayzo   = self._str_abs3D(yzo)

        xyz[0] = self.lat[0][0] + self.lat[1][0] + self.lat[2][0]
        xyz[1] = self.lat[0][1] + self.lat[1][1] + self.lat[2][1]
        xyz[2] = self.lat[0][2] + self.lat[1][2] + self.lat[2][2]
        axyz   = self._str_abs3D(xyz)

        for iat in xrange(self.num_atoms):
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        rat[0] = (self.cpos[iat][0] + i*self.lat[0][0] + j*self.lat[1][0] + k*self.lat[2][0])
                        rat[1] = (self.cpos[iat][1] + i*self.lat[0][1] + j*self.lat[1][1] + k*self.lat[2][1])
                        rat[2] = (self.cpos[iat][2] + i*self.lat[0][2] + j*self.lat[1][2] + k*self.lat[2][2])

                        projxoo = self._str_scl3D(rat, xoo) / axoo / axoo
                        projyoo = self._str_scl3D(rat, yoo) / ayoo / ayoo
                        projzoo = self._str_scl3D(rat, zoo) / azoo / azoo
                        projxyo = self._str_scl3D(rat, xyo) / axyo / axyo
                        projxzo = self._str_scl3D(rat, xzo) / axzo / axzo
                        projyzo = self._str_scl3D(rat, yzo) / ayzo / ayzo
                        projxyz = self._str_scl3D(rat, xyz) / axyz / axyz

                        if (\
                                (projxoo > -0.5 and projxoo <= 0.5) and\
                                (projyoo > -0.5 and projyoo <= 0.5) and\
                                (projzoo > -0.5 and projzoo <= 0.5) and\
                                (projxyo > -0.5 and projxyo <= 0.5) and\
                                (projxzo > -0.5 and projxzo <= 0.5) and\
                                (projyzo > -0.5 and projyzo <= 0.5) and\
                                (projxyz > -0.5 and projxyz <= 0.5)\
                           ):
                            self.cpos[iat][0] = rat[0]
                            self.cpos[iat][1] = rat[1]
                            self.cpos[iat][2] = rat[2]
                            i=10
                            j=10
                            k=10

    def setNames(self, in_names_vec):
        cnt               = -1
        in_names_vec_size = len(in_names_vec)

        for i in xrange(in_names_vec_size):
            if i < len(self.num_each_type):
                for j in xrange(self.num_each_type[i]):
                    cnt             += 1
                    self.names[cnt]  = in_names_vec[i]
                    if len(in_names_vec[i]) != 0:
                        self.names_were_given[cnt] = 1
                    else:
                        self.names_were_given[cnt] = 0


    def setAllAtomPos(self, in_names):
        self.names            = in_names
        self.names_were_given = np.array([0 for i in xrange(len(in_names))])

        for i in xrange(len(in_names)):
            if len(self.names[i]) > 0:
                self.names_were_given[i] = 1


    def setNamesWereGiven(self, in_names_were_given):
        self.names_were_given = in_names_were_given


    def setMom1(self, mom1_new):
        mom1_old = self.getMom1()
        dmom1    = uf.VVdiff(mom1_new, mom1_old)

        for iat in xrange(self.num_atoms):
            self.cpos[iat] = uf.SVprod(self.scale, self.cpos[iat])
            self.cpos[iat] = uf.VVsum(self.cpos[iat], dmom1)
            self.cpos[iat] = uf.SVprod(1.0 / self.scale, self.cpos[iat])

        self.calcAllDirectCoords()

    def setSD(self, in_sd):
        cnt  = -1
        size = len(in_sd)
        for i in xrange(size):
            if i < len(self.num_each_type):
                for j in xrange(self.num_each_type[i]):
                    cnt += 1
                    self.sd[cnt] = in_sd[i]

    def setISD(self, in_isd):
        self.isd = in_isd

    def addAtom(self, atom):
        at_type    = atom.getType()
        at_dpos    = atom.getDpos()
        at_cpos    = atom.getCpos()
        at_name    = atom.getName()

        at_name_given = 1
        if at_name == "NA":
            at_name_given = 0

        id = 0
        if at_type >= self.num_types:
            id                 = self.num_atoms
            self.num_each_type = np.append(self.num_each_type, [1])
            self.num_types    += 1
        else:
            for i in xrange(at_type):
                id += self.num_each_type[i]
            self.num_each_type[at_type] += 1

        np.insert(self.cpos, id, at_cpos)
        np.insert(self.dpos, id, at_dpos)
        np.insert(self.names, id, at_name)
        np.insert(self.names_were_given, id, at_name_given)
        self.num_atoms += 1


    def removeAtom(self, id):
        if id >= self.num_atoms:
            return

        np.delete(self.cpos, id, 0)
        np.delete(self.dpos, id, 0)
        np.delete(self.names, id, 0)
        np.delete(self.names_were_given, id, 0)

        self.num_each_type[self.types[id]] -= 1
        if self.num_each_type[self.types[id]] == 0:
            np.delete(self.num_each_type, types[id], 0)

        np.delete(self.types, id, 0)
        self.num_types = len(self.num_each_type)
        self.num_atoms -= 1


    def putInCell(self):
        TOL = 1e-15

        for i in xrange(self.num_atoms):
            for j in xrange(3):
                self.dpos[i][j] = self.dpos[i][j] - int(self.dpos[i][j])
                if self.dpos[i][j] < -TOL:
                    self.dpos[i][j] += 1
                if abs(self.dpos[i][j] < TOL):
                    self.dpos[i][j] = 0.0

        self.calcAllCartCoords()

    def atom_atom_distance(self, A, B):
        return np.sqrt( (A[0] - B[0]) * (A[0] - B[0]) +  (A[1] - B[1]) * (A[1] - B[1]) + (A[2] - B[2]) * (A[2] - B[2]) )

    def putInCompact(self):
        TOL = 1e-15

        for i in xrange(self.num_atoms):
            for j in xrange(3):
                self.dpos[i][j] = self.dpos[i][j] - int(self.dpos[i][j])
                if self.dpos[i][j] < -TOL:
                    self.dpos[i][j] += 1
                if abs(self.dpos[i][j] < TOL):
                    self.dpos[i][j] = 0.0

        for iat in xrange(len(self.dpos)):
            for ic in xrange(3):
                self.cpos[iat][ic] = 0.0
                for jc in xrange(3):
                    self.cpos[iat][ic] += self.dpos[iat][jc] * self.lat[jc][ic]

        adref1pos = np.array([0.0, 0.0, 0.0])
        adtstpos  = np.array([0.0, 0.0, 0.0])
        adtrgpos  = np.array([0.0, 0.0, 0.0])

        acref1pos = np.array([0.0, 0.0, 0.0])
        actstpos  = np.array([0.0, 0.0, 0.0])
        actrgpos  = np.array([0.0, 0.0, 0.0])

        for i in xrange(self.num_atoms):
            #scan all atoms to move except the first
            if i == 0:
                continue

            adtrgpos[0] = self.dpos[i][0]
            adtrgpos[1] = self.dpos[i][1]
            adtrgpos[2] = self.dpos[i][2]

            actrgpos[0] = self.cpos[i][0]
            actrgpos[1] = self.cpos[i][1]
            actrgpos[2] = self.cpos[i][2]

            min_bond    = 1.0e6

            #scan over all reference atoms
            for ii in xrange(i):
                adref1pos[0] = self.dpos[ii][0]
                adref1pos[1] = self.dpos[ii][1]
                adref1pos[2] = self.dpos[ii][2]

                #scan over all reference atoms
                for ic in xrange(3):
                    acref1pos[ic] = 0.0
                    for jc in xrange(3):
                        acref1pos[ic] += adref1pos[jc] * self.lat[jc][ic]

                #roll over first neighbour cells
                for i1 in [-1, 0, 1]:
                    for j1 in [-1, 0, 1]:
                        for k1 in [-1, 0, 1]:

                            adtstpos[0] = self.dpos[i][0] + i1
                            adtstpos[1] = self.dpos[i][1] + j1
                            adtstpos[2] = self.dpos[i][2] + k1

                            for ic in xrange(3):
                                actstpos[ic] = 0.0
                                for jc in xrange(3):
                                    actstpos[ic] += adtstpos[jc] * self.lat[jc][ic]
                            #test the bond distance
                            bond = self.atom_atom_distance(actstpos, acref1pos)
                            #if bond < min_bond, if it is OK DO IT
                            if ( bond < 1.03 * min_bond and abs(i - ii) < 10 ) or ( bond < 0.98 * min_bond ):
                                #update
                                min_bond    = bond

                                adtrgpos[0] = adtstpos[0]
                                adtrgpos[1] = adtstpos[1]
                                adtrgpos[2] = adtstpos[2]

                                actrgpos[0] = actstpos[0]
                                actrgpos[1] = actstpos[1]
                                actrgpos[2] = actstpos[2]

            self.dpos[i][0] = adtrgpos[0]
            self.dpos[i][1] = adtrgpos[1]
            self.dpos[i][2] = adtrgpos[2]

            self.cpos[i][0] = actrgpos[0]
            self.cpos[i][1] = actrgpos[1]
            self.cpos[i][2] = actrgpos[2]

        self.calcAllCartCoords()

    def shiftPos(self, shift, flag):
        #cartesian shift
        if flag == 0:
            for ia in xrange(self.num_atoms):
                self.cpos[ia] = uf.VVsum(shift, self.cpos[ia])
            self.calcAllDirectCoords()

        #direct coords shift
        if flag == -1:
            for ia in xrange(self.num_atoms):
                self.dpos[ia] = uf.VVsum(shift, self.dpos[ia])
            self.calcAllCartCoords()

    def RemoveCopies(self):
        tol    = 1e-3
        remove = False
        for ia1 in reversed(xrange(self.num_atoms)):
            for ia2 in xrange(ia1):
                if np.linalg.norm(uf.VVdiff(self.cpos[ia1], self.cpos[ia2])) < tol:
                    remove = True
            if remove:
                self.removeAtom(ia1)
            remove = False

    def rotate(self, rm):
        '''
        Rotation is done around the origin of the structure.
        Define origin q, a point p, and a rotation of p around
        q (R_q(p)).  Then R_q(p)=R_0(p-q)+q=R_0(p)-R_0(q)+q.
        We can evalutate the R_0 terms by simply mulitplying by
        a matrix.  Then we must add q to all the final cartesian
        positions.
        '''
        #Get R_0(p) for all cartesian positions.
        nlat     = np.transpose(self.lat)
        self.lat = np.transpose(rm * nlat)

        self.calcAllCartCoords()

        #Get R_0(q)
        r_orig   = rm * self.origin
        #Assign new cartesian positions
        for ia in xrange(self.num_atoms):
            self.cpos[ia] = uf.VVsum(uf.VVdiff(self.cpos[ia], r_orig), self.origin)
        #Get all the direct coords.
        self.calcAllDirectCoords()

    def setOrigin(self, in_origin):
        self.origin = in_origin

    def getTitle(self):
        return self.title

    def getScale(self):
        return self.scale

    def getLat(self):
        return self.lat

    def getScaledLat(self):
        tlat = np.matrix([[0.0, 0.0, 0.0] for i in xrange(3)])
        for ic in xrange(3):
            tlat[ic] = uf.SVprod(self.scale, self.lat[ic])
        return tlat

    def getNumEachType(self):
        return self.num_each_type

    def getTypes(self):
        return self.types

    def getType(self,id):
        return self.types[id]

    def getCoordFlag(self):
        return self.coord_flag

    def getCpos(self):
        return self.cpos

    def getDpos(self):
        return self.dpos

    def getNames(self):
        return self.names

    def getNumAtoms(self):
        return self.num_atoms

    def getAtomNameVec(self):
        return self.at_name_vec

    def getNamesWereGiven(self):
        return self.names_were_given

    def getAtomNumVec(self):
        return self.at_num_vec

    def getDispToAtomDir(self, in_dpos, at_num):
        #returns in_dpos-(dpos of atom at_num) in direct coordinates.  in_dpos should be in direct coords.
        disp = np.array([0.0, 0.0, 0.0])
        for ic in xrange(3):
            disp[ic] = in_dpos[ic] - self.dpos[at_num][ic]
        return disp

    def getDistToAtomDir(self, in_dpos, at_num):
        #returns norm(in_dpos-(pos of atom at_num)).  in_dpos should be in direct coords.
        disp = self.getDispToAtomDir(in_dpos, at_num)
        disp = uf.vecD2C(self.lat, disp)
        return self.scale * np.linalg.norm(disp)

    def getXrayScattFactor(in_name, in_lambda):
        #do not use lambda for now
        scatt_fact = 0.0
        for i in xrange(self.num_atoms):
            for j in xrange(len(self.at_name_vec)):
                if in_name == self.at_name_vec[j]:
                    scatt_fact = self.xray_scatt_vec[j]
        return scatt_fact

    def getAtomMass(self, at_num=None, name=None):
        if at_num is None and name is None:
            print "ERROR: GetAtomMass: Invalid Input parameters"
            return -1

        if at_num is not None:
            return self.at_mass_vec[at_num - 1]

        at_mass = -1
        for i in xrange(self.num_atoms):
            for j in xrange(len(self.at_name_vec)):
                if name == self.at_name_vec[j]:
                    at_mass = self.at_mass_vec[j]

        if at_mass < 0:
            print "ERROR: GetAtomMass: Cannot Find Atom with Name"

        return at_mass


    def getMom1(self):
        mom1 = np.array([0.0 for i in xrange(3)])
        for iat in xrange(self.num_atoms):
            mom1 = uf.VVsum(mom1, self.cpos[iat])

        if self.num_atoms > 0:
            mom1 = uf.SVprod(self.scale / (self.num_atoms * 1.0), mom1)
        else:
            mom1 = np.array([0.0 for i in xrange(3)])
        return mom1

    def getDispToAtomImageDir(self, in_dpos, at_num):
        '''
        returns {in_dpos}-{dpos of atom at_num} in direct coordinates, where the difference is taken
        between the nearest images (braces denote the fact that each set point is really a set of points
        in different image cells.  in_dpos should be in direct coords.
        '''
        disp = np.array([0.0, 0.0, 0.0])
        for ic in xrange(3):
            disp[ic] = in_dpos[ic] - self.dpos[at_num][ic]
            disp[ic] = disp[ic] - int(disp[ic])
            if disp[ic] > 0.5:
                disp[ic] -= 1
            if disp[ic] < -0.5:
                disp[ic] += 1
        return disp

    def getDistToAtomImageDir(self, in_dpos, at_num):
        '''
        returns {in_dpos}-{dpos of atom at_num} in direct coordinates, where the difference is taken
        between the nearest images (braces denote the fact that each set point is really a set of points
        in different image cells.  in_dpos should be in direct coords.
        '''
        disp = self.getDispToAtomImageDir(in_dpos, at_num)
        disp = uf.vecD2C(self.lat, disp)
        return self.scale * np.linalg.norm(disp)

    def getOrigin(self):
        return self.origin

    def calcAllCartCoords(self):
        for iat in xrange(len(self.dpos)):
            for ic in xrange(3):
                self.cpos[iat][ic] = 0.0
                for jc in xrange(3):
                    self.cpos[iat][ic] = self.cpos[iat][ic] + self.dpos[iat][jc] * self.lat[jc][ic]


    def calcAllDirectCoords(self):
        latT = self.lat.transpose()
        b    = self.cpos.transpose()
        if self.num_atoms > 0:
            latT, b = gaussJordan(latT, b)
            self.dpos = b.transpose()
        else:
            self.dpos = np.matrix([])
