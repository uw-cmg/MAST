import datetime
import time

def enum(**enums):
    return type('Enum',(),enums)

def enum(*sequential, **named):
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

def set2str(A,maxlen=20):

    if type(A) is not set:
        print 'INPUT IS NOT A SET.'
        return ""
    
    out = "{"
    lA = list(A)
    elements = ""
    nelem = len(lA)
    for j in range(nelem):
        if j < len(lA)-1:
            elements += '{}, '.format(lA[j])

    if nelem >0:
        elements += '{}'.format(lA[nelem-1])
        
    if len(elements) > maxlen-2:
        elements = elements[0:maxlen-5]+"..."
    out="{"+elements+"}"
    return out



# This part is for developing job and parser
def get_parent_child_sets(tokens):
    '''This function returns two sets as a tuple.
        Ex) (pset, cset) = get_parent_child_sets(tokens):
    '''
    parentset = set()
    childset = set()

    isparent = False
    ischild = False
    for item in tokens:
        if item.lower() == 'parent':
            isparent = True
            continue
        elif item.lower() == 'child':
            ischild = True
            isparent = False
            continue

        if isparent:
            parentset.add(item)
        if ischild:
            childset.add(item)
    return (parentset, childset)

def gettype(ingredientobj):
    stype=str(type(ingredientobj))
    stype = stype[1:-1]
    tokens = stype.replace("'","").split('.')
    return tokens[-1]


Jobstatus = enum('PreQ','InQ','Complete','Error')
JOB = Jobstatus

now = datetime.datetime.utcnow()
MAXSID = 10000
MAXJID = 100000
SCHEDULING_INTERVAL = 10
