# PickleManger example

from MAST.utility.picklemanager import PickleManager

pm = PickleManager()
A = 1
pm.save_variable(A)
A_loaded = pm.load_variable()
B = 'string'
class dummy:
    def __init__(self, id = None):
        if id is None:
            self.id = 1
        self.id = id
C =  dummy(123)

pm.save_variables(varlist=['A','B','C'], gdict=globals())

del A, B, C
# Before loading variables, the definition of class should be known.
# import class or find_class() of pickle should be run beforehand.
vardict = pm.load_variables()

# restore all variables from pickle to workspace.
pm.load_variables_to_ws(gdict=globals())



