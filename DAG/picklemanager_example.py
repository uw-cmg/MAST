# PickleManger example

from picklemanager import PickleManager
pm = PickleManager()
A = 1
B = 'string'
class dummy:
    def __init__(self, id = None):
        if id is None:
            self.id = 1
        self.id = id
C =  dummy(123)

pm.save_variables(['A','B','C'],globals())

del A, B, C
# Before loading variables, the definition of class should be known.
# import class or find_class() of pickle should be run beforehand.
vardict = pm.load_variables()

# restore all variables from pickle to workspace.
pm.load_variables_ws(gdict=globals())



