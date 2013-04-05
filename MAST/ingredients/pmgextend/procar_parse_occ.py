#!/usr/bin/env python

"""
Classes for reading/manipulating/writing VASP ouput files.
"""

from __future__ import division

__author__ = "Tam Mayeshiba"
__date__ = "April 4, 2013"

import os
import glob
import re
import math
import itertools
import warnings
import xml.sax.handler
import StringIO
from collections import defaultdict
import logging

import numpy as np

from pymatgen.util.coord_utils import get_points_in_sphere_pbc
from pymatgen.util.io_utils import zopen, clean_lines, micro_pyawk, \
    clean_json, reverse_readline
from pymatgen.core.structure import Structure
from pymatgen.core.composition import Composition
from pymatgen.electronic_structure.core import Spin, Orbital
from pymatgen.electronic_structure.dos import CompleteDos, Dos
from pymatgen.electronic_structure.bandstructure import BandStructure, \
    BandStructureSymmLine, get_reconstructed_band_structure
from pymatgen.core.lattice import Lattice
from pymatgen.io.vaspio.vasp_input import Incar, Kpoints, Poscar

logger = logging.getLogger(__name__)

import pdb



class ProcarOcc(object):
    """
    Object for reading a PROCAR file, only occupied bands
    """
    def __init__(self, filename):
        """
        Args:
            filename:
                Name of file containing PROCAR.
        """
        #create and return data object containing the information of a PROCAR
        self.name = ""
        self.data = dict()
        self.databands=""
        self.dataenergy=""
        self._read_file(filename)

    def get_d_occupation(self, atom_num):
        row = self.data[atom_num]
        return row[2]

    def _read_file(self, filename):
        with zopen(filename, "r") as f:
            lines = list(f.readlines()) #do not strip out comments
            #lines = list(clean_lines(f.readlines()))
        self.name = lines[0]
        kpointexpr = re.compile("^\s*k-point\s+(\d+).*weight = ([0-9\.]+)")
        expr = re.compile("^\s*([0-9]+)\s+")
        dataexpr = re.compile("[\.0-9]+")
        weight = 0
        myocc=0
        myband=0
        myenergy=0
        mystats=lines[1].split()
        numkpts = int(mystats[3])
        numbands = int(mystats[7])
        numions = int(mystats[11])
        self.databands = np.zeros([numions+1, numbands+1, 5],float)
        self.dataenergy = np.zeros([numions+1, numbands+1],float)
        #pdb.set_trace()
        for l in lines:
            if kpointexpr.match(l):
                m = kpointexpr.match(l)
                currentKpoint = int(m.group(1))
                weight = float(m.group(2))
                if currentKpoint == 1:
                    self.data = dict()
            if l.find("occ.") > -1:
                myocc=float(l.split()[7])
                myband=int(l.split()[1])
                myenergy=float(l.split()[4])
            if (myocc == 1) and expr.match(l):
                linedata = dataexpr.findall(l)
                linefloatdata = map(float, linedata)
                index = int(linefloatdata.pop(0))
                if index in self.data:
                    self.data[index] += np.array(linefloatdata) * weight
                else:
                    self.data[index] = np.array(linefloatdata) * weight
                self.databands[index][myband] += np.array(linefloatdata) * weight
                self.dataenergy[index][myband] += myenergy * weight

    def print_matrix(self):
        for key in self.data.iterkeys():
            print key,["{0:0.3f}".format(i) for i in self.data[key]]
        np.set_printoptions(precision=3, suppress=True, threshold=5000)
        print self.databands[9]
        dect = 0
        while dect < len(self.dataenergy[9]):
            print self.dataenergy[9][dect]
            dect = dect + 1


