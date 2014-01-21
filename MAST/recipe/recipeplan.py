############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
import os
import time
import logging
import inspect
import subprocess
from MAST.ingredients.chopingredient import ChopIngredient
from customlib.customchopingredient import CustomChopIngredient
#from MAST.ingredients.chopingredient import WriteIngredient
#from MAST.ingredients.chopingredient import IsReadyToRunIngredient
#from MAST.ingredients.chopingredient import RunIngredient
#from MAST.ingredients.chopingredient import IsCompleteIngredient
#from MAST.ingredients.chopingredient import UpdateChildrenIngredient
from MAST.utility import MASTFile
from MAST.utility import MASTError
from MAST.utility import loggerutils
class RecipePlan:
    """Contains the entire recipe plan.
        Attributes:
            self.name <str>: Recipe name
            self.ingredients <dict>: Dictionary of ingredients
                                     and their statuses
            self.update_methods <dict>: Dictionary of each
                  ingredient's update method, per child
                  [ingredient][child] = update method for child
            self.parents_to_check <dict>: Dictionary of parents
                  to check for completion before each ingredient
                  may be run
            self.run_methods <dict>: Dictionary of ingredient
                                     run methods
            self.write_methods <dict>: Dictionary of ingredient
                                       write methods
            self.ready_methods <dict>: Dictionary of ingredient
                                       ready methods
            self.complete_methods <dict>: Dictionary of ingred.
                                          complete methods
            self.ingred_input_options <dict>: Dictionary of 
                       ingredient input options
            self.status <str>: Recipe status
            self.working_directory <str>: Recipe working directory
            self.logger <logging logger>: Control-level logger
            self.recipe_logger <logging logger>: Recipe-level logger
    """
    def __init__(self, name, working_directory):
        self.name            = name
        self.ingredients     = dict()  #name: status
        self.update_methods   = dict()
        self.parents_to_check= dict()
        self.run_methods      = dict()
        self.write_methods    = dict()
        self.ready_methods    = dict()
        self.complete_methods = dict()
        self.ingred_input_options = dict()
        self.status="I"
        self.working_directory = working_directory
        self.logger = logging.getLogger('mastmon')
        self.logger = loggerutils.add_handler_for_control(self.logger)
        self.recipe_logger = logging.getLogger(self.working_directory)
        self.recipe_logger = loggerutils.add_handler_for_recipe(self.working_directory, self.recipe_logger)

    def do_ingredient_methods(self, iname, methodtype, childname=""):
        """Do the ingredient methods.
            Args:
                iname <str>: ingredient name
                methodtype <str>: method type (write, run, etc.)
        """
        if methodtype == 'mast_write_method':
            mdict = self.write_methods[iname]
        elif methodtype == 'mast_run_method':
            mdict = self.run_methods[iname]
        elif methodtype == 'mast_ready_method':
            mdict = self.ready_methods[iname]
        elif methodtype == 'mast_complete_method':
            mdict = self.complete_methods[iname]
        elif methodtype == 'mast_update_children_method':
            mdict = self.update_methods[iname][childname]
        else:
            raise MASTError(self.__class__.__name__,"Bad call to do_ingredient_methods with method type %s" % methodtype)
        allresults = list()
        self.logger.info("Do methods for %s" % methodtype)
        for methoditem in mdict.keys():
            minputs = list(mdict[methoditem])
            if methodtype == 'mast_update_children_method':
                minputs.append(childname)
            mresult = self.run_a_method(iname, methoditem, minputs)
            allresults.append(mresult)
        return allresults

    def run_a_method(self, iname, methodname, minputs):
        """Run a method. Evaluates which method to run
            using this order:
            A ChopIngredient method (no ".py")
            A CustomIngredient method (no ".py" and not found
                in ChopIngredient)
            A script-defined method (some "def xxxx" in customlib/*.py)
            Args:
                iname <str>: ingredient name
                methodname <str>: method name
                minputs <list>: method inputs
            Returns:
                Results of the method, whatever it returns.
        """
        len_inputs = len(minputs)
        #Is it a ChopIngredient method?
        handlerlist = inspect.getmembers(ChopIngredient)
        for htry in handlerlist:
            hname = htry[0]
            if methodname == hname: #found it. Return mresult.
                my_ing = ChopIngredient(name = self.ingred_input_options[iname]['name'],
                    program_keys=self.ingred_input_options[iname]['program_keys'],
                    structure=self.ingred_input_options[iname]['structure'])
                if len_inputs == 0:
                    mresult = getattr(ChopIngredient, methodname)(my_ing)
                elif len_inputs == 1:
                    mresult = getattr(ChopIngredient, methodname)(my_ing, minputs[0])
                elif len_inputs == 2:
                    mresult = getattr(ChopIngredient, methodname)(my_ing, minputs[0], minputs[1])
                elif len_inputs == 3:
                    mresult = getattr(ChopIngredient, methodname)(my_ing, minputs[0], minputs[1], minputs[2])
                elif len_inputs == 4:
                    mresult = getattr(ChopIngredient, methodname)(my_ing, minputs[0], minputs[1], minputs[2], minputs[3])
                else:
                    raise MASTError(self.__class__.__name__, "Function %s for ChopIngredient requires too many inputs (> 4)." % methodname)
                self.recipe_logger.info("Results for method name %s: %s" % (methodname, mresult))
                return mresult
        #Is it a CustomIngredient method?
        handlerlist = inspect.getmembers(CustomChopIngredient)
        for htry in handlerlist:
            hname = htry[0]
            if methodname == hname: #found it. Return mresult.
                my_ing = CustomChopIngredient(name = self.ingred_input_options[iname]['name'],
                    program_keys=self.ingred_input_options[iname]['program_keys'],
                    structure=self.ingred_input_options[iname]['structure'])
                if len_inputs == 0:
                    mresult = getattr(CustomChopIngredient, methodname)(my_ing)
                elif len_inputs == 1:
                    mresult = getattr(CustomChopIngredient, methodname)(my_ing, minputs[0])
                elif len_inputs == 2:
                    mresult = getattr(CustomChopIngredient, methodname)(my_ing, minputs[0], minputs[1])
                elif len_inputs == 3:
                    mresult = getattr(CustomChopIngredient, methodname)(my_ing, minputs[0], minputs[1], minputs[2])
                elif len_inputs == 4:
                    mresult = getattr(CustomChopIngredient, methodname)(my_ing, minputs[0], minputs[1], minputs[2], minputs[3])
                else:
                    raise MASTError(self.__class__.__name__, "Function %s for ChopIngredient requires too many inputs (> 4)." % methodname)
                self.recipe_logger.info("Results for method name %s: %s" % (methodname, mresult))
                return mresult
        #Is it some script in customlib?
        mpath = os.path.join(os.getenv("MAST_INSTALL_PATH"),"customlib",methodname)
        if os.path.isfile(mpath):
            subinputs = list()
            subinputs.append("python")
            subinputs.append(mpath)
            subinputs.extend(minputs)
            myproc = subprocess.Popen(subinputs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            myproc.wait()
            (mresult, merr) = myproc.communicate()
            self.recipe_logger.info("Results for method name %s: %s %s" % (methodname, mresult, merr))
            return mresult
        else:
            raise MASTError(self.__class__.__name__,"Could not find method or .py file %s" % methodname)
        return None




    def write_ingredient(self, iname):
        """Write the ingredient files according to the 
            correct method
        """
        return self.do_ingredient_methods(iname, "mast_write_method")

    def complete_ingredient(self, iname):
        """Check if an ingredient is complete
        """
        cresults = self.do_ingredient_methods(iname, "mast_complete_method")
        iscomplete = None
        for cresult in cresults:
            if (cresult == None) or (cresult == ""):
                pass
            elif cresult == False:
                return False
            elif cresult == True:
                iscomplete = True
            else:
                pass
        self.recipe_logger.info("Completeness evaluation: %s (from %s)" % (iscomplete, cresults))
        return iscomplete

    def ready_ingredient(self, iname):
        """Check if an ingredient is ready
        """
        return self.do_ingredient_methods(iname, "mast_ready_method")


    def run_ingredient(self, iname):
        """Run ingredient
        """
        return self.do_ingredient_methods(iname, "mast_run_method")

    def update_children(self, iname):
        """Update the children of an ingredient
        """
        upd_results=list()
        for childname in self.update_methods[iname].keys():
            upd_results.append(self.do_ingredient_methods(iname, "mast_update_children_method", childname))
        return upd_results

    def fast_forward_check_complete(self):
        """Check if runs are complete."""
        for iname in self.ingredients.keys():
            if not (self.ingredients[iname] in ["C", "E", "skip"]):
                if self.complete_ingredient(iname):
                    self.ingredients[iname] = "C"
                    self.recipe_logger.info("Status of %s changed to %s" % (iname, "C"))
                    self.update_children(iname)
        return


    def check_if_have_parents(self):
        """Check if runs at "Initialized" status have parents 
            and switch them to "Wait" if so; otherwise,
            switch them to "Stage"
        """
        for iname in self.ingredients.keys():
            if self.ingredients[iname] == "I":
                ptc = list(self.parents_to_check[iname])
                plen = len(ptc)
                if plen > 0:
                    self.ingredients[iname] = "W"
                    self.recipe_logger.info("Status of %s changed to %s" % (iname, "W"))
                else:
                    self.ingredients[iname] = "S"
                    self.recipe_logger.info("Status of %s changed to %s" % (iname, "S"))


    def check_if_ready_to_proceed_are_complete(self):
        """Check if ready-to-proceed ingredients are complete
        """
        for iname in self.ingredients.keys():
            if self.ingredients[iname] == "P":
                if self.status_change_recommended(iname):
                    pass
                else:
                    if self.complete_ingredient(iname):
                        self.ingredients[iname] = "C"
                        self.update_children(iname)
        return
    
    def status_change_recommended(self, iname):
        """Check if a status change is recommended for the ingredient,
            as listed in the ingredient folder/change_status.txt.
            Args:
                iname <str>: ingredient name
            Returns:
                True if a status change was recommended, and 
                    changes the status of the ingredient in self.ingredients.
                False otherwise
        """
        statuspath = os.path.join(self.working_directory, iname, "change_status.txt")
        if not os.path.isfile(statuspath):
            return False
        statusfile = MASTFile(statuspath)
        newdata=list()
        changed=False
        for sline in statusfile.data: #status:recommend:timestamp
            if not "status_changed" in sline:
                newstatus = sline.split(":")[0]
                self.ingredients[iname]=newstatus
                newline = sline + ":status_changed:" + time.asctime() + "\n"
                self.recipe_logger.info("Status of %s changed to %s" % (iname, newstatus))
                changed=True
                newdata.append(newline)
            else:
                newdata.append(sline)
        statusfile.data=list(newdata)
        statusfile.to_file(statuspath)
        return changed

    def check_if_parents_are_complete(self):
        """Check if parents of waiting ingredients are
            complete.
        """
        for iname in self.ingredients.keys():
            if self.ingredients[iname] == "W":
                okay=0
                ptc = list(self.parents_to_check[iname])
                plen = len(ptc)
                for parent in ptc:
                    if self.ingredients[parent] == "C":
                        self.recipe_logger.info("Status of %s changed to %s" % (parent, "C"))
                        okay = okay + 1
                if okay == plen:
                    self.ingredients[iname] = "S"
                    self.recipe_logger.info("Status of %s changed to %s" % (iname, "S"))

    def run_staged_ingredients(self):
        """Run staged ingredients.
        """
        for iname in self.ingredients.keys():
            if self.ingredients[iname] == "S":
                if self.complete_ingredient(iname):
                    self.ingredients[iname] = "C"
                    self.recipe_logger.info("Status of %s changed to %s" % (iname, "C"))
                    self.update_children(iname)
                else:
                    if not (self.ready_ingredient(iname)):
                        self.write_ingredient(iname)
                    if self.ready_ingredient(iname):
                        self.run_ingredient(iname)
                        self.ingredients[iname] = "P"
                        self.recipe_logger.info("Status of %s changed to %s" % (iname, "P"))

    def check_recipe_status(self, verbose=1):
        """Check ingredient statuses, and get recipe status
            I = Initialized
            W = Waiting on parents
            S = Staged
            P = ready to Proceed
            C = Complete
        """
        #self.fast_forward_check_complete()
        self.check_if_have_parents()
        self.check_if_ready_to_proceed_are_complete()
        self.check_if_parents_are_complete()
        self.run_staged_ingredients()
        self.print_status(verbose)
    def print_status(self, verbose=1):
        """Print status and set the recipe status.
        """

        total = len(self.ingredients.keys())
        totcomp=0
        totwait=0
        totproceed=0
        totinit=0
        totstage=0
        toterr=0
        totskip=0
        ilist = self.ingredients.keys()
        ilist.sort()
        statusfile = MASTFile()
        statusfile.data.append("#Statuses:\n")
        statusfile.data.append("#I = Initial state\n")
        statusfile.data.append("#W = Waiting on parents\n")
        statusfile.data.append("#S = Staged for run prep (parents complete)\n")
        statusfile.data.append("#P = ready to Proceed to queue; look for submission\n")
        statusfile.data.append("#C = Complete\n")
        statusfile.data.append("#E = Error\n")
        statusfile.data.append("#skip = Skip (from user)\n")


        if verbose == 1:
            namestring = "Recipe name: %s" % self.name
            self.recipe_logger.info(namestring)
            self.recipe_logger.info(time.asctime())
            self.logger.info(namestring)
            self.logger.info(time.asctime())
        for iname in ilist:
            if verbose == 1:
                ingstring = "%30s : %4s" % (iname, self.ingredients[iname])
                self.recipe_logger.info(ingstring)
                self.logger.info(ingstring)
            statusfile.data.append("%30s : %4s\n" % (iname, self.ingredients[iname]))
            if self.ingredients[iname] == "C":
                totcomp = totcomp + 1
            elif self.ingredients[iname] == "P":
                totproceed = totproceed + 1
            elif self.ingredients[iname] == "I":
                totinit = totinit + 1
            elif self.ingredients[iname] == "W":
                totwait = totwait + 1
            elif self.ingredients[iname] == "S":
                totstage = totstage + 1
            elif self.ingredients[iname] == "E":
                toterr = toterr + 1
            elif self.ingredients[iname] == "skip":
                totskip = totskip + 1
        headerstring = "%8s %8s %8s %8s %8s %8s %8s= %8s" % ("INIT","WAITING","STAGED","PROCEED","COMPLETE", "ERROR", "USERSKIP", "TOTAL")
        self.recipe_logger.info(headerstring)
        self.logger.info(headerstring)
        valuestring = "%8i %8i %8i %8i %8i %8i %8i= %8i" % (totinit, totwait, totstage, totproceed, totcomp, toterr, totskip, total)
        self.recipe_logger.info(valuestring)
        self.logger.info(valuestring)
        if totcomp == total:
            self.status = "C"
        else:
            self.status = "R"
        #print "Recipe status: %s" % self.status
        statusfile.to_file(os.path.join(self.working_directory,"status.txt"))
        if toterr > 0:
            self.logger.error("ATTENTION: Recipe at %s has one or more ingredient errors. Please check the MAST_ERROR file and or mast_recipe.log file, in the recipe directory. Change status 'E' in status.txt to 'W' when the error has been resolved." % self.working_directory)

    def add_ingredient(self, ingredient_name, ingredient):
        """Used to add an ingredient_object corresponding to an ingredient name
        """
        raise MASTError(self.__class__.__name__, "This function is obsolete.")
        self.ingredients[ingredient_name] = ingredient

    def get_ingredient(self, ingredient_name):
        """Used to get an ingredient_object corresponding to an ingredient name
        """
        raise MASTError(self.__class__.__name__, "This function is obsolete.")
        return self.ingredients.get(ingredient_name)

    #def add_parent(self, ingredient_name, parent_name):
    #    """Used to build the dependency dict using the parent child relationship
    #    """
    #    self.dependency_dict.setdefault(ingredient_name, list()).append(parent_name)

    def ingredient_iterator(self):
        """Iterates through the ingredients dict and returns the ingredients one by one
        """
        raise MASTError(self.__class__.__name__, "This function is obsolete.")
        for ingredient_name, ingredient_obj in self.ingredients.iteritems():
            yield ingredient_name, ingredient_obj

    def __iter__(self):
        """Iterates through the ingredients dict and returns the ingredients one by one
        """
        raise MASTError(self.__class__.__name__, "This function is obsolete.")
        for ingredient_name, ingredient_obj in self.ingredients.iteritems():
            yield ingredient_obj
    def __repr__(self):
        """Print information."""
        rlines=""
        rlines=rlines + "Recipe name: %s\n" % self.name
        rlines=rlines + "Ingredients: \n"
        ikeys = self.ingredients.keys()
        ikeys.sort()
        for ikey in ikeys:
            rlines=rlines + "    %s\n" % ikey
            rlines=rlines + "        Parents: \n" 
            for myparent in self.parents_to_check[ikey]:
                rlines=rlines + "            %s\n" % myparent
            rlines=rlines + "        Write:    %s\n" % self.write_methods[ikey]
            rlines=rlines + "        Ready:    %s\n" % self.ready_methods[ikey]
            rlines=rlines + "        Run:      %s\n" % self.run_methods[ikey]
            rlines=rlines + "        Complete: %s\n" % self.complete_methods[ikey]
            rlines=rlines + "        Children: \n"
            for htukey in self.update_methods[ikey]:
                rlines=rlines + "            %s:%s\n" % (htukey, self.update_methods[ikey][htukey])
        self.print_status()
        return rlines
    def get_statuses_from_file(self):
        """Get status of each ingredient from a status.txt
            file in the recipe directory.
        """
        statpath = os.path.join(self.working_directory,'status.txt')
        if not os.path.isfile(statpath):
            raise MASTError(self.__class__.__name__, "Could not get status file in %s" % statpath)
        statfile = MASTFile(statpath)
        for statline in statfile.data:
            if statline[0] == "#":
                pass
            else:
                statsplit = statline.strip().split(':')
                oneingred = statsplit[0].strip()
                onestatus = statsplit[1].strip()
                if oneingred in self.ingredients.keys():
                    self.ingredients[oneingred] = onestatus
                else:
                    raise MASTError(self.__class__.__name__, "Ingredient %s is not in the original recipe's ingredients list." % oneingred)

            
