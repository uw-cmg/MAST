def ofile(name):
    l=open(name,'a')
    return l

def restart_output(Optimizer):
    outputoptions = ['output', 'summary','Genealogyfile', 'tenergyfile', 'debugfile','fpfile', 'fpminfile']
    outfiles = {'files':None, 'ifiles':None, 'output':None, 'summary':None,
        'Genealogyfile':None, 'tenergyfile':None, 'debugfile':None,
        'fpfile':None, 'fpminfile':None}
    for one in outputoptions:
        outpar = eval('Optimizer.{0}'.format(one))
        if outpar:
            outfiles[one] = eval("open({0},'a')".format(repr(outpar)))
    outfiles['files'] = map(ofile,Optimizer.files)
    if Optimizer.indiv_defect_write:
        outfiles['ifiles'] = map(ofile,Optimizer.ifiles)
    return outfiles
