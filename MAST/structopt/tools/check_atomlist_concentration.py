import copy
import logging

def check_atomlist_concentration(atomlist,natoms, loggername=None):
    """Checks format of atomlist with natoms to ensure number and
    structure of atoms is correct"""
    flag = False
    checkal = sum([abs(c) for ind,c,m,u in atomlist])
    natoms = int(natoms)
    if loggername:
        #logger = initialize_logger(loggername)
        logger = logging.getLogger(loggername)
        logger.info('Check_atomlist_concentration recieved atomlist={0} and natoms={1}'.format(
            atomlist, natoms))
    if checkal < 1:
        if loggername:
            logger.critical('Number of atoms in atomlist cannot be less than 1.')
            logger.debug('Number of atoms is : {0}'.format(checkal))
        raise RuntimeError('Error: In Atom List, Concentration <1')
        # There cannot be less than one atom in the system
        catmlist = None
    elif checkal == 1:
        if natoms != 1:
            if loggername:
                logger.info('Assuming atom list concentration given as fraction of total atoms')
            natmlist = [[atomlist[i][0],int(round(natoms*atomlist[i][1])),atomlist[i][2],atomlist[i][3]] for i in range(len(atomlist))]
            flag = True
        else:
            catmlist = copy.copy(atomlist)
    elif checkal == 100:
        if natms != 100:
            if loggername:
                logger.info('Assuming atom list concentrations give as percent of total atoms')
            natmlist=[(atomlist[i][0],int(round(natoms*atomlist[i][1]/100)),atomlist[i][2],atomlist[i][3]) for i in range(len(atomlist))]
            flag = True
        else:
            catmlist = copy.copy(atomlist)
    else:
        catmlist = copy.copy(atomlist)
    if flag:
        atms=sum([c for ind,c,m,u in natmlist])
        extraatms=natoms-atms
        while extraatms > 0:
            random.choice(natmlist)[1]+=1
            atms=sum([c for ind,c,m,u in natmlist])
            extraatms=natoms-atms
        while extraatms < 0:
            random.choice(natmlist)[1]-=1
            atms=sum([c for ind,c,m,u in natmlist])
            extraatms=natoms-atms
        catmlist = [tuple(one) for one in natmlist]
    return catmlist