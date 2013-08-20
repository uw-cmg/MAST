#!/usr/bin/env python
from MAST.utility.picklemanager import PickleManager as pm
import sys

mypm = pm(sys.argv[1])
vardict = mypm.load_variable()
print vardict
