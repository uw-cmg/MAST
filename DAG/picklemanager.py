import pickle

class PickleManager:
    def __init__(self, filename=None):
        self._protocol = 2 # binary
        if filename is None:
            self.fn = 'mast.pkl'
        else:
            self.fn = filename

    def save_variables(self, varlist,gdict):
        '''gdict = globals()'''
        fh = open(self.fn,'wb')
        pickle.dump(varlist,fh)
        for key in varlist:
            pickle.dump(gdict[key],fh, self._protocol)
        fh.close()

    def load_variables(self, filename=None):
        '''Default filename is initialized when object is initialized.'''
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

    def load_variables_to_ws(self, filename=None, gdict=globals()):
        ''' Load variables to workspace from a pickle file.
        Usage 1 : load_variables_to_ws(self, gdict = globals())
        Usage 2 : load_variables_to_ws(self, filename = 'savefile.pkl',gdict = globals())
        '''
        if filename is None:
            fn = self.fn
        else:
            fn = filename
        fh = open(fn,'rb')
        varlist = pickle.load(fh)
        for key in varlist:
            if key in gdict:
                print '%s is overwritten.' % key
            gdict[key] = pickle.load(fh)



