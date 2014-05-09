##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
raise NotImplementedError("Comment out this line if you are sure you want to run this file! You may want to delete initialize.rst before running make latexpdf or make html")
import os
#  cmd = "sphinx-apidoc --separate --force -o ./source %s" % os.environ['MAST_INSTALL_PATH']
cmd = "sphinx-apidoc --separate -o ./source %s" % os.environ['MAST_INSTALL_PATH']
os.system(cmd)
