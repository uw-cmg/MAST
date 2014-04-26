##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import time
from MAST.utility import dirutil
from MAST.utility import MASTError

class MASTFile:
    """Controls MAST file IO and manipulation"""
    #data = [] #TTM 11/16/11 remove from here to avoid confusing reinit.

    #TTM 10/31/11 allow to initialize without a file path
    def __init__(self, file_path=str()):
        self.data = list() # TTM 11/16/11 add here
        self.p = dict() # TTM 12/22/11 set in file class instead #dictionary holding a mapping of param names to param values
        if file_path != "":
            if not os.path.isfile(file_path):
                raise MASTError(self.__class__.__name__,
                "No file at %s; cannot create MASTFile" % file_path)
            self.from_file(file_path)
        else:
            pass

    def from_file(self, file_path):
        """Reads data from a file and stores it in a list"""
        self.data = []
        #TTM+2 10/7/11 add error checking in case of no file
        if not os.path.isfile(file_path):
            raise MASTError(self.__class__.__name__,
                "No such file at " + file_path)
        readf = open(file_path,'rb')
        lines = readf.readlines()
        for line in lines:
            self.data.append(line)
        readf.close()    
    
    #TTM 12/6/11 create for multiple species potcars
    def from_file_append(self, file_path=str()):
        """Reads data from a file and appends it to the list"""
        if not os.path.isfile(file_path):
            return
        readf = open(file_path,'rb')
        lines = readf.readlines()
        for line in lines:
            self.data.append(line)
        readf.close()    
    
    def to_file(self,file_path):
        """Writes data to a file (overwrites existing file)."""
        #TTM+2 10/7/11 add error checking in case of no path
        if (file_path == "") or (file_path == None):
            raise MASTError(self.__class__.__name__,
                "Refusing to copy to empty path.")
        if self.data == []:
            raise MASTError(self.__class__.__name__,
                "Empty file not copied to " + file_path)
        dirutil.lock_directory(os.path.dirname(file_path))
        writef = open(file_path,'wb')
        for line in self.data:
            writef.write(line)
        writef.close()
        dirutil.unlock_directory(os.path.dirname(file_path))

    #TTM 11/10/11 created to solve archiving problems
    def to_unique_file(self, parent_path="", try_name="", suffix="", max=10):
        """Writes to a unique file (?)"""
        if (parent_path == "") or (parent_path == None):
            return
        if (try_name == "") or (try_name == None):
            try_name = str(time.time())
        basicname = os.path.join(parent_path, try_name)
        tryme = basicname + suffix
        ct = 1
        while os.path.isfile(tryme):
            tryme = basicname + str(ct) + suffix
            ct = ct + 1
            if ct > max: #TTM 11/23/11 increase max threshold
                print 'Over', max, 'of the same file. Erroring out.'
                return
        self.to_file(tryme)
        return tryme #TTM 20121109 return the file name

    def to_stdout(self):
        """Writes to stdout"""
        for line in self.data:
            print line

    #TTM+block 10/5/11 created (Pgetline)
    def get_line_number(self, line_number):
        """Grabs what line is at number line_number"""
        matidx = line_number - 1 # matrix indices start at 0
        if matidx < 0:
            return None
        if len(self.data) > matidx:
            return self.data[line_number - 1]
        else:
            return None

    #TTM+block 10/5/11 created (Pgetline)
    def get_line_match(self, string_to_match):
        """Searches for a string, and returns it once found"""
        for line in self.data:
            if string_to_match in line:
                return line
        return None

    #TTM+block 11/7/11 created to return line number
    def get_line_match_number(self, string_to_match):
        linect = len(self.data)
        ct = 0
        myline=""
        while ct < linect:
            myline = self.data[ct]
            if string_to_match in myline:
                return ct + 1 #return line numbers starting from 1
            ct = ct + 1

    #TTM+block 10/5/11 created
    def get_last_line_match(self, string_to_match):
        lastmatch=""
        for line in self.data:
            if string_to_match in line:
                lastmatch = line
        if lastmatch == "":
            return None
        return lastmatch

    def get_last_x_lines_line_match(self, string_to_match, x_lines):
        lastmatch=""
        x_int = int(x_lines)
        if x_int == 0:
            return self.get_last_line_match(string_to_match)
        searchme = self.data[-1*x_int:]
        for line in searchme:
            if string_to_match in line:
                lastmatch = line
        if lastmatch == "":
            return None
        return lastmatch

    #TTM+block 11/8/11 separated from get_segment_from_last_line_match
    def get_segment(self, string_to_chop, start_string="", end_string=""):
        start_found=""
        startpos=""
        endpos=""
        if start_string == "":
            startpos = 0 # start from beginning
        else:
            start_found = string_to_chop.find(start_string)
            if start_found == -1:
                return None # not found
            startpos = start_found + len(start_string) # start after the string
        if end_string == "":
            endpos = len(string_to_chop) # go to end
        else:
            endpos = string_to_chop.find(end_string)
            if endpos == -1:
                return None #not found
        return string_to_chop[startpos:endpos].strip() # strip off white space


    #TTM+block 10/5/11 created (Pgetseg)
    def get_segment_from_last_line_match(self, string_to_match, start_string="", end_string="", match_last_line=True):
        found=""
        if match_last_line:
            found = self.get_last_line_match(string_to_match)
        else:
            found = self.get_line_match(string_to_match)
        if not found:
            return None
        myseg = self.get_segment(found, start_string, end_string)
        return myseg

    #TTM 11/7/11 from Pmodfilenew.py
    def modify_file_by_line_number(self, lineno="", mode="", param=""):
        """Modifies a file at a line number"""
        #mode: D - Delete, R - Replace, I - insert
        #param: if mode is:
        #       D - ignored; leave blank
        #       R - string with which to replace whole line at lineno
        #       I - string to insert below lineno as a whole line
        try:
            lineno = int(lineno)
        except (ValueError, TypeError):
            print 'Invalid line number.'
            raise RuntimeError('Line number does not exist')

        try:
            mode = str.upper(mode)
        except (ValueError, TypeError):
            print 'Invalid mode.'
            raise RuntimeError('Selected mode does not exist')

        if mode not in ['D','R','I']:
            print 'Only D-delete, R-replace, I-insert supported.'
            return

        dataidx = lineno - 1 #sync line number to matrix elements
        origlen = len(self.data)
        before=[]
        after=[]

        if (dataidx > origlen - 1):
            print 'File does not contain', lineno, 'lines.'
            return

        if (dataidx < 0):
            print 'Negative line number cannot be found.'
            return
        #
        if (dataidx == 0):
            before=[]
        else:
            before = self.data[0:dataidx]

        if (dataidx == origlen - 1):
            after=[]
        else:
            after = self.data[dataidx+1:]

        if (mode == 'D'):
            pass # no special action
        elif (mode == 'R'):
            before.append(param)# put new line where old one went
        elif (mode == 'I'):
            before.append(self.data[dataidx]) # put line back in
            before.append(param)

        self.data = before + after
        return
    
    #TTM+block 11/8/11 created
    #TTM 12/22/11 moving to file class from INCAR class
    def file_to_dictionary(self, from_path):
        self.from_file(from_path)
        line = ""
        linesplit = list()
        for line in self.data:
            if line == "":
                pass #empty space
            if not ('=' in line):
                pass # system comment
            if line[0] == "#":
                pass #comments
            linesplit = line.split('=',1) # split in 1 equals sign
            if len(linesplit) > 1:
                self.p[linesplit[0].strip()] = linesplit[1].strip()

    #TTM 052112 added
    def copy_data_to(self, other_file):
        other_file.data = list()
        for line in self.data:
            other_file.data.append(line)
        return
