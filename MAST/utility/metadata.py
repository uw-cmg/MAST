##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Glen Jenness
# Last updated:
#       2016-07-20 Tam Mayeshiba
#       2013-07-01
##############################################################
import os

from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility.mastfile import MASTFile

class Metadata(MASTFile):
    """Class to handle the metadata file
        Fields are composed of a keyword and a value, separated by a '='.
    """
    def __init__(self, metafile=""):
        self.mdpath=metafile
        if not (os.path.isfile(self.mdpath)):
            MASTFile.__init__(self, "")
        else:
            MASTFile.__init__(self, self.mdpath)
        return

    def write_data(self, keyword, data, option=0):
        """Writes a keyword and its associated data to the metafile"""
        [existline, existdata] = self.search_data(keyword)
        if existline == None:
            linetowrite = "%s = %s\n" % (keyword, data)
            self.data.append(linetowrite)
        else:
            if (option == 0): #append to previous entry
                linetowrite = "%s = %s; %s\n" % (keyword, existdata, data)
            elif (option == 1): #overwrite previous entry
                linetowrite = "%s = %s\n" % (keyword, data)
            self.modify_file_by_line_number(existline, "R", linetowrite)
        self.to_file(self.mdpath)
        return

    def search_data(self, keyword):
        """Searches the file for a keyword, and if found returns the line number
            and data for that keyword.
            Returns:
                [line_number, data]
                line_number <int>: line number found
                data <str>: data found
        """
        line_number = None
        data = None
        
        searchstr = keyword.lower() + " = "
        line_number = self.get_line_match_number(searchstr) #TTM use mastfile
        if not (line_number == None): #TTM match found
            dataidx = line_number - 1 #TTM adjust line number to data index
            linetext = self.data[dataidx]
            data = linetext.split(" = ")[1].strip()

        return line_number, data

    def read_data(self, keyword):
        """Searches the metadata file for a specific keyword and returns the
            data.
        """
        line_number, data = self.search_data(keyword)
        return data
 
    def clear_data(self, keyword):
        """Removes the specified data and keyword from the file.
        """
        line_number, data = self.search_data(keyword)
        if line_number == None: #TTM not found
            return
        self.modify_file_by_line_number(line_number, "D")
        self.to_file(self.mdpath)
        return

    def clear_file(self):
        """Empties the metafile"""
        self.data = list()
        self.to_file(self.mdpath)
        return

    def __repr__(self):
        self.to_stdout()
        return
