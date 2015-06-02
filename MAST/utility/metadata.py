##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Glen Jenness
# Last updated: 2013-07-01
##############################################################
import os

from MAST.utility import MASTObj
from MAST.utility import MASTError

ALLOWED_KEYS = {
                'metafile' : (str, 'metadata.txt', 'Metadata file name')
               }


class Metadata(MASTObj):
    """Class to handle the metadata file
        Fields are composed of a keyword and a value, separated by a '='.
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)

    def write_data(self, keyword, data, option=0):
        """Writes a keyword and its associated data to the metafile"""
        with open(self.keywords['metafile'], 'a') as metafile:
            # First check to see if the keyword already exists in the metadata file
            if None in self.search_data(keyword):
                metafile.write('%s = %s\n' % (keyword, data))
            else:
                if (option == 0):
                    entry = self.read_data(keyword)
                    entry += '; %s' % data
                    self.clear_data(keyword)
                    metafile.write('%s = %s\n' % (keyword, entry))
                elif (option == 1):
                    self.clear_data(keyword)
                    metafile.write('%s = %s\n' % (keyword, data))

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

        with open(self.keywords['metafile'], 'r') as metafile:
            for n, line in enumerate(metafile):
                if (keyword.lower() + ' = ') in line.lower():
                    line_number = n
                    data = line.split(' = ')[1].strip()
                    break

        return line_number, data

    def read_data(self, keyword):
        """Searches the metadata file for a specific keyword and returns the
            data.
        """
        line_number, data = self.search_data(keyword)
        return data
 
    def clear_data(self, keyword):
        """Removes the specified data and keyword from the file.
            This is inefficient, as what we do is basically reading in the file
            then iterating through the file, do if..else check, clear the file,
            then re-write it.

            But it works!
        """
        line_number, data = self.search_data(keyword)

        data = list()
        with open(self.keywords['metafile'], 'r') as metafile:
            for n, line in enumerate(metafile):
                if n == line_number:
                    pass
                else:
                    data.append(line)

        self.clear_file()
        for line in data:
            keyword, data = line.split(' = ')
            self.write_data(keyword, data.strip())

    def clear_file(self):
        """Empties the metafile"""
        open(self.keywords['metafile'], 'w').close()

    def __repr__(self):
        return file(self.keywords['metafile']).read()
