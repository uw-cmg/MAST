import random

def swap_vacancy(indiv, Optimizer):
    """Move function to perform mutation to swap the vacancy position with another atom position
    Inputs:
        indiv = Individual class object to be altered
        Optimizer = Optimizer class object with needed parameters
    Outputs:
        indiv = Altered Individual class object
    """
    if 'MU' in Optimizer.debug:
        debug = True
    else:
        debug = False
    Optimizer.output.write('Vacancy Swap Mutation performed on individual\n')
    Optimizer.output.write('Index = '+repr(indiv.index)+'\n')
    vacancies = indiv.vacancies.copy()
    if len(vacancies) == 0:
        Optimizer.output.write('Vacancy Swap Failed. Individual has no identified vacancies to exchange\n')
        natomsswap = 0
    elif len(vacancies) == 1:
        natomsswap = 1
    else:
        natomsswap = random.randint(1,len(vacancies))
    str=[]
    for one in range(natomsswap):
        vacant = random.choice(vacancies)
        if Optimizer.alloy:
            solid=indiv[0].copy()
            solid.extend(indiv.bulki)
            select = random.choice([one.index for one in solid if one.symbol==vacant.symbol])
        else:
            select = random.choice([one.index for one in indiv[0] if one.symbol==vacant.symbol])
        if select < len(indiv[0]):
            indiv[0][select].position, vacancies[vacant.index].position = vacant.position, indiv[0][select].position
            str.append('C')
        else:
            indiv.bulki[select-len(indiv[0])].position, vacancies[vacant.index].position = vacant.position, indiv.bulki[select-len(indiv[0])].position
            str.append('B')
    indiv.vacancies = vacancies
    Optimizer.output.write('Exchanged the position of '+repr(natomsswap)+' vacancies\n')
    cex = len([one for one in str if one=='C'])
    bex = len([one for one in str if one=='B'])
    Optimizer.output.write('Cluster exchanges: '+repr(cex)+'\n')
    Optimizer.output.write('Bulk exchanges: '+repr(bex)+'\n')
    muttype='VS'+repr(natomsswap)+'c'+repr(cex)+'b'+repr(bex)
    if indiv.energy==0:
        indiv.history_index=indiv.history_index+'m'+muttype
    else:
        indiv.history_index=repr(indiv.index)+'m'+muttype
    return indiv