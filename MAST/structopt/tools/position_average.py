def position_average(atoms, mode='median', trimp = 0.05):
    """Function to calculate the trimmed mean or median
    of the positions of atoms
    Input:
        atoms = ASE Atoms object
        mode = method to use either median or trimmedmean
            default value is median
        trimp = percentage of high and lows to leave out
            default value is 5% == 0.05
    Output:
        [x,y,z] = position of atoms center
    """
    positions = atoms.get_positions()
    xs = []
    ys = []
    zs = []
    for x,y,z in positions:
        xs.append(x)
        ys.append(y)
        zs.append(z)
    xs.sort()
    ys.sort()
    zs.sort()
    if mode=='median':
        xc = xs[len(xs)/2]
        yc = ys[len(ys)/2]
        zc = zs[len(zs)/2]
    else:
        ndel = int(round(trimp*len(xs)))
        xc = sum(xs[ndel:len(xs)-ndel])/(len(xs)-2*ndel)
        yc = sum(ys[ndel:len(ys)-ndel])/(len(ys)-2*ndel)
        zc = sum(zs[ndel:len(zs)-ndel])/(len(zs)-2*ndel)
    return [xc,yc,zc]