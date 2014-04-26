##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
#  cmd = "sphinx-apidoc --separate --force -o ./source %s" % os.environ['MAST_INSTALL_PATH']
cmd = "sphinx-apidoc --separate -o ./source %s" % os.environ['MAST_INSTALL_PATH']
os.system(cmd)
