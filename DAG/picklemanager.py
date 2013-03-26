import pickle
# later pickle to cPickle

class PickleManager:

    def __init__(self, filename=None):
        if filename is None:
            self.fn = 'mast.pkl'
        else:
            self.fn = filename
    "gdict = globals()"
    def save_variables(self, varlist,gdict):
        fh = open(self.fn,'wb')
        pickle.dump(varlist,fh)
        for key in varlist:
            pickle.dump(gdict[key],fh)
        fh.close()

    "Default filename is initialized when object is initialized."
    def load_variables(self, filename=None):
        vardict = {}
        if filename is None:
            fn = self.fn
        else:
            fn = filename
        fh = open(fn,'rb')
        varlist = pickle.load(fh)
        for key in varlist:
            vardict[key] = pickle.load(fh)
        
        return vardict
            
        


