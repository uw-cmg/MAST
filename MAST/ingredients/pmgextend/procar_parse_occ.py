#!/usr/bin/env python
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# This procar code should be replaced with pymatgen functions. 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################

"""
Classes for reading/manipulating/writing VASP ouput files.
"""

from __future__ import division

__author__ = "Tam Mayeshiba from vasp_output by Shyue Ping Ong"
__date__ = "April 8, 2012"

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

logger = logging.getLogger(__name__)




class ProcarOcc(object):
    """
    Object for reading a PROCAR file
    """
    def __init__(self, filename):
        """
        Args:
            filename:
                Name of file containing PROCAR.
        """
        #create and return data object containing the information of a PROCAR
        self.name = ""
        self.data = {}
        self.headers = None
        myocc=0
        with zopen(filename, "r") as f:
            #lines = list(clean_lines(f.readlines())) #TTM do not clean lines
            lines = list(f.readlines())
            self.name = lines[0]
            kpointexpr = re.compile("^\s*k-point\s+(\d+).*weight = ([0-9\.]+)")
            ionexpr = re.compile("^ion.*")
            expr = re.compile("^\s*([0-9]+)\s+")
            dataexpr = re.compile("[\.0-9]+")
            weight = 0

            for l in lines:
                if l.find("occ.")>-1:
                    myocc=float(l.split()[7])
                if kpointexpr.match(l):
                    m = kpointexpr.match(l)
                    currentKpoint = int(m.group(1))
                    weight = float(m.group(2))
                    if currentKpoint == 1:
                        self.data = dict()
                elif ionexpr.match(l) and self.headers is None:
                    self.headers = l.split()
                    self.headers.pop(0)
                elif expr.match(l):
                    linedata = dataexpr.findall(l)
                    linefloatdata = map(float, linedata)
                    index = int(linefloatdata.pop(0))
                    if index in self.data:
                        if myocc==1:
                            self.data[index] += np.array(linefloatdata) * weight
                    else:
                        if myocc==1:
                            self.data[index] = np.array(linefloatdata) * weight
                        else:
                            self.data[index]=0.0

    def get_d_occupation(self, atom_index):
        """
        .. deprecated:: v2.6.4

            Use get_occpuation instead.

        Returns the d occupation of a particular atom.

        Args:
            atom_index:
                Index of atom in PROCAR.

        Returns:
            d-occupation of atom at atom_index.
        """
        warnings.warn("get_d_occupation has been deprecated. Use "
                      "get_occupation instead.", DeprecationWarning)
        return self.get_occupation(atom_index, 'd')

    def get_occupation(self, atom_index, orbital):
        """
        Returns the occupation for a particular orbital of a particular atom.

        Args:
            atom_num:
                Index of atom in the PROCAR
            orbital:
                A string representing an orbital. If it is a single
                character, e.g., s, p, d or f, the sum of all s-type,
                p-type, d-type or f-type orbitals occupations are returned
                respectively. If it is a specific orbital, e.g., px, dxy,
                etc., only the occupation of that orbital is returned.

        Returns:
            Sum occupation of orbital of atom.
        """
        row = self.data[atom_index]
        total = 0
        found = False
        for orb, data in zip(self.headers, row):
            if orb.startswith(orbital):
                found = True
                total += data
        if not found:
            raise ValueError("Invalid orbital {}".format(orbital))
        return total


