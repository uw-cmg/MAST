##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-06-13 by Zhewen Song
##############################################################
import os, math

from MAST.utility import MASTObj
from MAST.utility import InputOptions
from MAST.utility import MASTError
from MAST.utility import Metadata
from MAST.utility import dirutil

ALLOWED_KEYS = {\
                 'templateFile'    : (list, None, 'template file name'),\
                 'inputOptions'    : (InputOptions, None, 'input options parsed using input parser'),\
                 'personalRecipe'  : (str, None, 'personalized recipe file'),\
                 'working_directory' : (str, None, 'Working directory'),
               }

class RecipeTemplateParser(MASTObj):
    """Class for parsing the template recipe file in the input file and
        creating a "personalized" recipe file where any <sys> and <N> from
        the recipe template have been filled in according to information
        from the input file.
        Attributes:
            self.input_options <InputOptions>: Input options from input file.
            self.template_file <str>: Path to the recipe template file
            self.personal_recipe <str>: Path to the personalized recipe.
            self.ingredient_list <list of str>: List of ingredients mentioned
                                                in the recipe template file.
            self.chunks <list of list>: List of chunks
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.input_options   = self.keywords['inputOptions']
        self.template_file   = self.keywords['templateFile']
        self.personal_recipe = self.keywords['personalRecipe']
        self.ingredient_list = list()

        self.metafile = Metadata(metafile='%s/metadata.txt' % self.keywords['working_directory'])
        self.chunks = list()
    def parse(self):
        """ Parses the template recipe file and creates
            the personalized recipe file
        """
        if len(self.template_file) == 0:
            raise MASTError(self.__class__.__name__, "Recipe contents not provided!")

        if self.input_options is None:
            raise MASTError(self.__class__.__name__, "Input Options not provided!")
        if self.personal_recipe is None:
            raise MASTError(self.__class__.__name__, "Personal recipe file not provided!")
        
        #fetch required paramaters
        #f_ptr           = open(self.template_file, "r")
        recipe_contents = list(self.template_file)
        #print recipe_contents
        o_ptr           = open(self.personal_recipe, "a")
        system_name     = self.input_options.get_item("mast", "system_name", "sys")
        n_defects       = self.input_options.get_item("defects", "num_defects", 0)
        d_defects       = self.input_options.get_item("defects","defects")
        n_images        = self.input_options.get_item("neb", "images", 0)
        d_neblines      = self.input_options.get_item("neb", "neblines", {})
        recipe_name     = None

#        print system_name, self.input_options.get_item('mast', 'system_name')
        chunkcount=0
        mychunk=list()
        modchunk=False
        for line in recipe_contents: #f_ptr.readlines():
            duplicate = line
            if ((len(duplicate) - len(duplicate.lstrip(' '))) % 4 != 0):
                raise MASTError("parsers/recipetemplateparser", "Recipe at %s contains incorrect number of whitespace chars at the beginning of the line! Please convert all indentations to the appropriate number of groups of four spaces." % line)
            if '\t' in line:
                raise MASTError("parsers/recipetemplateparser","The tab character exists in recipe template %s. Please convert all tabs to the appropriate number of groups of four spaces." % self.template_file)
            if '{begin}' in line:
                self.chunks.append(list(mychunk))
                mychunk=list()
            elif '{end}' in line:
                pass
            else:
                mychunk.append(line)
        if len(mychunk) > 0:
            self.chunks.append(list(mychunk))
        #for chunk in self.chunks:
        #    print chunk

        input_options_keys = self.input_options.get_sections()
        key = 'personal_recipe'
        expandedlist=list()
        if key in input_options_keys:
            return expandedlist
        else:
            for chunk in self.chunks:
                expanded=self.parse_chunk(chunk)
                expandedlist.extend(expanded)
            o_ptr.write("$personal_recipe\n")
            o_ptr.writelines(expandedlist)
            o_ptr.write("$end\n")
            #f_ptr.close()
            o_ptr.close()
            return expandedlist
        #self.chunks[chunkcount]=dict()
        #self.chunks[chunkcount]['modify'] = modchunk
        #    #line             = line.strip()
        #    #line             = line.lower()
        #    processing_lines = []
        #    #shortcut straight copy + 6
        #    processing_lines.append(line)
        #    output_str = "\n".join(processing_lines)
        #    o_ptr.write("%s\n" % output_str)

        #f_ptr.close()
        #o_ptr.close()
#        print 'in RecipeParser.parse():', list(set(self.ingredient_list))
        #return recipe_name
    def parse_chunk(self, chunk):
        """Parse a chunk of lines.
            Args:
                chunk <list>: List of lines
            Returns:
                expandedchunk <list>: List of lines
        """
        origchunk = list(chunk)
        expandedchunk = list()
        needsdefects=0
        needscharges=0
        needsphonons=0
        needsnebs=0
        needsscaling=0
        for line in chunk:
            if "<N>" in line:
                needsdefects=1
            if "<B-E>" in line:
                needsnebs=1
            if "<P>" in line:
                needsphonons=1
            if "<Q>" in line:
                needscharges=1
            if "<B>" in line:
                needsnebs=1
            if "<E>" in line:
                needsnebs=1
            if "<S>" in line:
                needsscaling=1
        d_scaling       = self.input_options.get_item("scaling")
        d_defects       = self.input_options.get_item("defects","defects")
        d_nebs          = self.input_options.get_item("neb","nebs")
        
        if needsscaling == 1:
            scalingsize = d_scaling.keys()
            scalingsize.sort()
        else: scalingsize = ['1x1x1']

        for size in scalingsize:
            if needsdefects == 1:
                mydefects=d_defects.keys()
                mydefects.sort()
                for defectname in mydefects:
                    for charge in d_defects[defectname]['charge']:
                        if charge < 0:
                            mycharge = 'q=n' + str(int(math.fabs(charge)))
                        else:
                            mycharge = 'q=p' + str(int(charge))
                        if needsphonons == 1:
                            if len(d_defects[defectname]['phonon'].keys()) > 0:
                                phononkeys = d_defects[defectname]['phonon'].keys()
                                phononkeys.sort()
                                for phonon in phononkeys:
                                    for line in origchunk:
                                        newline = line.replace("<N>", defectname)
                                        if needscharges == 1:
                                            newline = newline.replace("<Q>", mycharge)
                                        if needsscaling == 1:
                                            newline = newline.replace("<S>",size)
                                        newline = newline.replace("<P>", phonon)
                                        expandedchunk.append(newline)
                        else:
                            for line in origchunk:
                                newline = line.replace("<N>", defectname)
                                if needscharges == 1:
                                    newline = newline.replace("<Q>", mycharge)
                                if needsscaling == 1:
                                    newline = newline.replace("<S>",size)
                                expandedchunk.append(newline)
            elif needsnebs == 1:
                nebkeys = d_nebs.keys()
                nebkeys.sort()
                for neblabel in nebkeys:
                    defbegin = neblabel.split('-')[0]
                    defend = neblabel.split('-')[1]
                    chargebegin = d_defects[defbegin]['charge']
                    chargeend = d_defects[defend]['charge']
                    chargeboth = set(chargebegin) & set(chargeend)
                    for charge in chargeboth:
                        if charge < 0:
                            mycharge = 'q=n' + str(int(math.fabs(charge)))
                        else:
                            mycharge = 'q=p' + str(int(charge))
                        if needsphonons == 1:
                            if len(d_nebs[neblabel]['phonon'].keys()) > 0:
                                phononkeys = d_nebs[neblabel]['phonon'].keys()
                                phononkeys.sort()
                                for phonon in phononkeys:
                                    for line in origchunk:
                                        newline = line.replace("<B>", defbegin)
                                        newline = newline.replace("<E>", defend)
                                        newline = newline.replace("<B-E>", neblabel)
                                        if needscharges == 1:
                                            newline = newline.replace("<Q>", mycharge)
                                        if needsscaling == 1:
                                            newline = newline.replace("<S>",size)
                                        newline = newline.replace("<P>", phonon)
                                        expandedchunk.append(newline)
                        else:
                            for line in origchunk:
                                newline = line.replace("<B>", defbegin)
                                newline = newline.replace("<E>", defend)
                                newline = newline.replace("<B-E>", neblabel)
                                if needscharges == 1:
                                    newline = newline.replace("<Q>", mycharge)
                                if needsscaling == 1:
                                    newline = newline.replace("<S>",size)                            
                                expandedchunk.append(newline)
            elif needsscaling==1:
                for line in origchunk:
                    newline = line.replace("<S>",size)
                    expandedchunk.append(newline)
            else: expandedchunk = list(origchunk)
        return expandedchunk
        #origchunk = list(expandedchunk)
        #expandedchunk=list()
        #for defectname in self.d_defects:
        #    for line in origchunk:
        #        newline = line.replace("<N>", defectname)
        #        expandedchunk.append(line)


    def old_parsing(self):
        raise MASTError(self.__class__.__name__, "This function is obsolete.") 
        for line in linestr:
            #validate the input line
            if not line or line.startswith('#'):
                continue

            #collect recipe name
            line2 = line.split()
            if (line2[0].lower() == 'recipe'):
                recipe_name = line2[1] 

            #collect ingredents
            line2 = line.split()
            if (line2[0].lower() == 'ingredient'):
                self.ingredient_list.append(line2[2])


            #replace line params<N>
            #step 1 : replace <sys> with system_name
            #step 2 : replace <n-n> with appropriate hop combinations
            #step 3 : replace <img-n> with appropriate image numbers
            #step 4 : replace <n> with appropriate defect numbers
            processing_lines.append(line)
            #step 1
            processing_lines = self.process_system_name(processing_lines, system_name)
            #step 4
            processing_lines = self.process_defects(processing_lines, n_defects, d_defects)
            #step 2
            processing_lines = self.process_hop_combinations(processing_lines, d_neblines)
            #step 3
            processing_lines = self.process_images(processing_lines, n_images)

            self.make_metadata_entries(processing_lines)

            self.process_phononlines(processing_lines)
            #dump the processed lines to file
            output_str = "\n".join(processing_lines)
            o_ptr.write("%s\n" % output_str)

        f_ptr.close()
        o_ptr.close()
#        print 'in RecipeParser.parse():', list(set(self.ingredient_list))
        return recipe_name

    def process_system_name(self, processing_lines, system_name):
        """replace <sys> with the system name from the input options
        """
        raise MASTError(self.__class__.__name__, "This function is obsolete.") 
        for index in xrange(len(processing_lines)):
            processing_lines[index] = processing_lines[index].replace('<sys>', system_name)
        return processing_lines

    def process_hop_combinations(self, processing_lines, d_neblines):
        """replace <n-n> with neb labels which are keys of the
           neblines dict of the input options.

            Args:
                processing_lines <list of str>: recipe lines to process.
                d_neblines <dict of str>: dictionary of NEB lines.
        """
        raise MASTError(self.__class__.__name__, "This function is obsolete.") 
        new_lines = []
        eval_lines = []
        if not d_neblines:
            return processing_lines
        line=""
        for line in processing_lines:
            if "<n-n>" in line:
                for neblabel in d_neblines.keys():
                    n_line = line.replace('<n-n>', neblabel)
                    eval_lines.append(n_line)
                    if 'ingredient' in n_line and not '<' in n_line: 
                        keyword = n_line.split()[1]
                        data = 'neblabel: %s' % neblabel
                        self.metafile.write_data(keyword, data)
            else:
                eval_lines.append(line)
        line=""
        for line in eval_lines:
            #print "TTM DEBUG line: ", line
            if not 'neb' in line.split('_'):
                #print "TTM DEBUG line safe"
                new_lines.append(line)
            else:
                evalsplit = line.split('child')
                if len(evalsplit) == 1:
                    new_lines.append(line)
                else:
                    childsplit=evalsplit[1].split('_')
                    parsplit=evalsplit[0].split('_')
                    okay=1
                    for neblabel in d_neblines.keys():
                        if okay == 1 and (neblabel in childsplit):
                            parlabels = neblabel.split('-')
                            #print "TTM DEBUG: parlabels: ",parlabels
                            #print "TTM DEBUG: parsplit: ",parsplit
                            if not parlabels[0] in parsplit and not (parlabels[1] in parsplit):
                                if not neblabel in parsplit: #image static
                                    okay=0
                    if okay == 0:
                        pass
                    else:
                        new_lines.append(line)
        return new_lines

    def process_images(self, processing_lines, n_images):
        """replace <img-n> with the equivalent number of 
           lines based on the number of images found in the
           input_options
        """
        raise MASTError(self.__class__.__name__, "This function is obsolete.") 
        new_lines = []
        if not n_images:
            return processing_lines
        for line in processing_lines:
            if '<img-n>' in line:
                for index in xrange(n_images):
                    n_line = line.replace('<img-n>', str(index+1).zfill(2))
                    if 'ingredient' in n_line and not '<' in n_line: 
                        keyword = n_line.split()[1]
                        data = 'name: %s' % keyword
                        self.metafile.write_data(keyword, data)
                    new_lines.append(n_line)
            else:
                 new_lines.append(line)
        return new_lines

    def process_defects(self, processing_lines, n_defects, d_defects):
        """replace <N> with the equivalent number of lines
           based on the number of defects given in the
           input options

           Args:
            processing_lines <list of str>: recipe lines to proceess
            n_defects <int>: number of defected systems
            d_defects <dict>: dictionary of defects, including labels and 
                                positions.
        """
        raise MASTError(self.__class__.__name__, "This function is obsolete.") 
        #import inspect
        #print 'GRJ DEBUG: %s.%s' % (self.__class__.__name__, inspect.stack()[0][3])
        #print d_defects

        #print 'GRJ DEBUG: parse_defects() working_directory =', self.keywords['working_directory']
        new_lines = list()

        if not n_defects:
            return processing_lines

        for line in processing_lines:
            #print 'GRJ DEBUG: line =', line
            if ('<n>' in line) or ('<q>' in line):
                for defect_key in d_defects.keys():
                    #print 'GRJ DEBUG: defect_key =', defect_key
                    #defect_label = defect_key.split('_')[1] #defect_1, etc.
                    def_line = line.replace("<n>", defect_key)

                    charge_list = d_defects[defect_key]['charge']
                    #print 'GRJ DEBUG: charge_list =', charge_list
                    #print 'GRJ DEBUG: def_line before charges:', def_line
                    for charge in charge_list:
                        if (charge < 0):
                            clabel = 'q=n' + str(abs(charge))
                        else:
                            clabel = 'q=p' + str(charge)
                        #print 'GRJ DEBUG: clabel =', clabel
                        new_def_line = def_line.replace('<q>', clabel)
                        new_lines.append(new_def_line)
                        #new_lines.append(def_line.replace('<q>', clabel))

                        if 'ingredient' in def_line:
                            keyword = new_def_line.split()[1]
                            data = 'defect_label: %s, charge: %i' % (defect_key, charge)
                            #print 'GRJ DEBUG: def_line, keyword, charge', new_def_line, keyword, charge
                            self.metafile.write_data(keyword, data)
            else:
                new_lines.append(line)

        return new_lines

    def process_phononlines(self, processing_lines):
        """add phonon information to the metadata. Does not change line info.
        """
        raise MASTError(self.__class__.__name__, "This function is obsolete.") 
        for line in processing_lines:
            if 'ingredient' in line and 'phonon_' in line:
                nameval = line.split()[1]
                [dataline,dataval]=self.metafile.search_data(nameval)
                okay=0
                if not (dataval == None):
                    datapcs = dataval.split(',')
                    for datapc in datapcs:
                        dlabel = datapc.split(":")[0].strip()
                        dval = datapc.split(":")[1].strip()
                        if dlabel == 'neblabel' or dlabel == 'defect_label':
                            data = 'phononlabel: %s' % dval
                            self.metafile.write_data(nameval, data)
                            okay=1
                            break
                if okay==0:
                    if 'perfect' in line:
                        data = 'phononlabel: perfect'
                        self.metafile.write_data(nameval, data)
                        okay=1
                    else:
                        data = 'phononlabel: %s' % nameval
                        self.metafile.write_data(nameval, data)
                        okay=1
        return

    def make_metadata_entries(self, processing_lines):
        """Add metadata entry for all ingredients. 
            Does not change line information.
        """
        raise MASTError(self.__class__.__name__, "This function is obsolete.") 
        for line in processing_lines:
            if 'ingredient' in line:
                nameval = line.split()[1]
                data = 'name: %s' % nameval
                self.metafile.write_data(nameval, data)
        return 
    def get_unique_ingredients(self):
        """fetches the ingredients names"""
        raise MASTError(self.__class__.__name__, "This function is obsolete.") 
        return list(set(self.ingredient_list))
