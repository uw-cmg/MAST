import os
cmd = "sphinx-apidoc -f -o ./source %s" % os.environ['MAST_INSTALL_PATH']
os.system(cmd)
