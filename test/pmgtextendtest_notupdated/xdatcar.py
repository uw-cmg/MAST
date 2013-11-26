from MAST.ingredients.pmgextend import vasp_extensions
import os
testdir=os.path.join(os.getenv("MAST_INSTALL_PATH"),"test/pmgextend_test")
myxdat = vasp_extensions.read_my_xdatcar(testdir,"XDATCAR_vasp5211")
myxdatold = vasp_extensions.read_my_xdatcar(testdir,"XDATCAR_vasp522")
vasp_extensions.write_my_xdatcar(testdir,myxdat,"written_5211")
vasp_extensions.write_my_xdatcar(testdir,myxdatold,"written_522")
