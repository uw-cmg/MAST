##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import logging
import os
import time
from MAST.utility import dirutil

def get_mast_logger(loggername=""):
    """Get the MAST logger. Add a handler if none exists.
        Args:
            loggername <str>: Logger name
    """
    if loggername == "":
        logger = logging.getLogger()
    else:
        logger = logging.getLogger(loggername)
    logger = add_mast_handler(logger)
    return logger

def add_mast_handler(logger, hname="mast.log"):
    """Add a handler and format string
        Args:
            logger <logging.Logger>: Logger for which to add the handler
            hname <str>: handler name (file name)
    """
    logger.setLevel(logging.INFO)
    if os.getenv("MAST_DEBUG") == None:
        pass
    else:
        logger.setLevel(logging.DEBUG) 
    formatstr ='%(levelname)8s: %(name)s: %(asctime)s: %(message)s'
    formatter = logging.Formatter(formatstr)
    monitorhandler = logging.FileHandler(filename=os.path.join(dirutil.get_mast_control_path(), hname))
    monitorhandler.setFormatter(formatter)
    #controlfilter = ControlFilter()
    #controlhandler.addFilter(controlfilter)
    if not getattr(logger, 'has_monitor_handler', None):
        logger.addHandler(monitorhandler)
        logger.has_monitor_handler = True
    return logger


def validate_recipe_name(rname):
    """Validate the recipe name.
        Use $MAST_CONTROL if recipe name is not found in $MAST_SCRATCH.
        Args:
            rname <str>: Recipe name candidate
        Returns:
            rnamevalid <str>: Validated recipe name
    """
    #check existence of rname:
    mast_scratch = dirutil.get_mast_scratch_path()
    rdirs = os.listdir(mast_scratch)
    rnamevalid = dirutil.get_mast_control_path()
    for rdir in rdirs:
        if rdir in rname:
            rnamevalid = os.path.join(mast_scratch, rdir)
            return rnamevalid
    return rnamevalid

def config_for_recipe_and_control(rname):
    rnamevalid = validate_recipe_name(rname)
    formatstr ='%(asctime)s : %(module)10s: %(levelname)8s : %(message)s'
    longformatter = logging.Formatter(formatstr)
    formatstrshort ='%(levelname)8s : %(message)s'
    shortformatter = logging.Formatter(formatstrshort)
    recipehandler = logging.FileHandler(filename="%s/mast_recipe.log" % rnamevalid)
    recipehandler.setFormatter(longformatter)
    controlhandler = logging.FileHandler(filename="%s/mast.log" % dirutil.get_mast_control_path())
    controlhandler.setFormatter(shortformatter)
    logger=logging.getLogger('')
    if not getattr(logger, 'has_recipe_handler', None):
        logger.addHandler(recipehandler)
        logger.has_recipe_handler = True
    if not getattr(logger, 'has_control_handler', None):
        logger.addHandler(controlhandler)
        logger.has_control_handler = True
    logger.setLevel(logging.INFO)
    return logger

def add_handler_for_recipe(rname, logger):
    """Add a handler for the recipe.
        Args:
            rname <str>: Recipe full path
            logger <logging.Logger>: Logger
    """
    #print "Handler added for rname %s at %s" % (rname, time.asctime())
    logger.setLevel(logging.INFO)
    rnamevalid = validate_recipe_name(rname)
    formatstr ='%(asctime)s : %(module)10s: %(name)20s: %(levelname)8s : %(message)s'
    formatter = logging.Formatter(formatstr)
    recipehandler = logging.FileHandler(filename="%s/mast_recipe.log" % rnamevalid)
    recipehandler.setFormatter(formatter)
    rnameshort = os.path.basename(rnamevalid)
    recipefilter = RecipeFilter(rnameshort)
    recipehandler.addFilter(recipefilter)
    if not getattr(logger, 'has_recipe_handlers', None):
        logger.has_recipe_handlers = list()
    if not rnameshort in logger.has_recipe_handlers:
        logger.addHandler(recipehandler)
        logger.has_recipe_handlers.append(rnameshort)
    return logger

def add_handler_for_control(logger):
    """Add a handler for control
        Args:
            logger <logging.Logger>: Logger
    """
    logger.setLevel(logging.INFO)
    #formatstr ='%(levelname)8s : %(name)10s: %(message)s'
    formatstr ='%(levelname)8s : %(message)s'
    formatter = logging.Formatter(formatstr)
    controlhandler = logging.FileHandler(filename="%s/mast.log" % dirutil.get_mast_control_path())
    controlhandler.setFormatter(formatter)
    controlfilter = ControlFilter()
    controlhandler.addFilter(controlfilter)
    if not getattr(logger, 'has_control_handler', None):
        logger.addHandler(controlhandler)
        logger.has_control_handler = True
    return logger

class ControlFilter(logging.Filter):
    def filter(self, record):
        return True
        if getattr(record, 'name', None) in ['mast','mastmon']:
            return True
        else:
            return False

class RecipeFilter(logging.Filter):
    def __init__(self, rnameshort):
        """Recipe filter
            Args:
                rnameshort <str>: Recipe name (not full path)
        """
        logging.Filter.__init__(self)
        self.rname = rnameshort

    def filter(self, record):
        if self.rname in getattr(record, 'name', None):
            return True
        else:
            return False
        

def initialize_logger(filename="default.log", formatstr=""):
    logger     = logging.getLogger(filename)
    if not getattr(logger, 'handler_set', None):
        handler    = logging.FileHandler(filename)
        if formatstr == "":
            #formatter  = logging.Formatter('%(asctime)s : %(module)15s:%(lineno)4d> : %(levelname)8s : %(message)s')
            formatter  = logging.Formatter('%(asctime)s : %(module)10s: %(levelname)8s : %(message)s')
        else:
            formatter = logging.Formatter(formatstr)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        logger.handler_set = True
    return logger

def initialize_short_logger(filename="default.log"):
    format = "%(levelname)8s: %(message)s"
    #format = "%(asctime)s %(levelname)8s : %(message)s"
    return initialize_logger(filename, format)
