from MAST.ingredients.performneb import PerformNEB

class PerformNEBLowMesh(PerformNEB):
    def __init__(self, **kwargs):
        PerformNEB.__init__(self, **kwargs)
