def remove_duplicates(oldlist, tolerance):
    """Function to remove the duplicates in a list.
    ***A duplicate is defined as being less than a given tolerance
    different from the other numbers in the list***
    Inputs:
        oldlist = List of floats/integers possibly containing duplicates
        tolerance = maximum difference between items for them to be considered duplicates
    Outputs:
        newlist = List without duplicates
        newlistindices = List of indices of items in newlist corresponding to positions in old list
    """
    newlist = []
    dups = []
    newlistindices = []
    for i in range(len(oldlist)):
        diffs = []
        #Identify the difference between given item and other items in list
        for j in range(len(oldlist)):
            diffs.append(abs(oldlist[i]- oldlist[j]))
        #Check if any of the differences are less than the tolerance
        flaggedlist = []
        for k in range(len(diffs)):
            if diffs[k] <= tolerance:
                if k !=i:
                    flaggedlist.append(k)
        if len(flaggedlist) != 0:
            dups.append([i,flaggedlist])
        else:
            newlistindices.append(i)
    selects = []
    for one,others in dups:
        if len(selects)==0:
            selects.append(one)
        else:
            flag=True
            for two in others:
                if two in selects:
                    flag=False
            if flag:
                selects.append(one)
    newlistindices.extend(selects)
    newlistindices.sort()
    for one in newlistindices:
        newlist.append(oldlist[one])
    
    return newlist, newlistindices