#-------------------------------------------------------------------------------
# Name:        Atom
# Purpose:
#
# Author:      Kumaresh visakan
#
# Created:     11/01/2013
#-------------------------------------------------------------------------------

import numpy as np

class Atom:
    #constructor initialized with default values
    def __init__(self):
        self.name         = "NA"
        self.num          = -1
        self.unit_cell    = np.array([0, 0, 0])
        self.cpos         = np.array([-999.0, -999.0, -999.0])
        self.dpos         = np.array([-999.0, -999.0, -999.0])
        self.corigin      = np.array([0.0, 0.0, 0.0])
        self.type         = -1

    #setter methods for member variables
    def setName(self, name):
        self.name = name

    def setNum(self, num):
        self.num = num

    def setUnitCell(self, unit_cell):
        if isinstance(unit_cell, np.array) and len(unit_cell) == 3:
            self.unit_cell = unit_cell

    def setCpos(self, cpos):
        if isinstance(cpos, np.array) and len(cpos) == 3:
            self.cpos = cpos

    def setDpos(self, dpos):
        if isinstance(dpos, np.array) and len(dpos) == 3:
            self.dpos = dpos

    def setCOrigin(self, corigin):
        if isinstance(corigin, np.array) and len(corigin) == 3:
            self.corigin = corigin

    def setType(self, in_type):
        self.type = in_type


    #getter methods for member variables
    def getName(self):
        return self.name

    def getNum(self):
        return self.num

    def getUnitCell(self):
        return self.unit_cell

    def getCpos(self):
        return self.cpos

    def getDpos(self):
        return self.dpos

    def getCOrigin(self):
        return self.corigin

    def getType(self):
        return self.type

    #get displacement from origin
    def getCDispFromOrigin(self):
        diff = np.array([0, 0, 0])
        for index in xrange(3):
            diff[index] = self.cpos[index] - corigin[index]
        return diff

    #get distance from origin
    def getDistFromOrigin(self):
        return np.linalg.norm(self.getCDispFromOrigin())


