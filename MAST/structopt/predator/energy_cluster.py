from MAST.structopt.tools import get_best
from MAST.structopt.predator import mutation_dups
import random

def energy_cluster(pop, Optimizer):
    """Use a k-means clustering approach to identify individuals for new population
    """
    #Get coordinates for each individual
    pts = [ind.energy for ind in pop]
    #Identify spatial bounds
    mins = min(pts)
    maxs = max(pts)
    STR=''
    attemptcount=10
    passflag=False
    while attemptcount > 0:
        #Select random initial centroid locations
        cents = []
        while len(cents) < Optimizer.nindiv:
            a = random.choice(pts)
            if a not in cents:
                cents.append(a)
        #Assign individual to closest centroid
        clustlist = [[] for i in range(Optimizer.nindiv)]
        for one in pts:
            ds = []
            for i in range(len(cents)):
                d = abs(one - cents[i])
                ds.append([d,i])
            ds.sort()
            loc = ds[0]
            clustlist[loc[1]].append(one)
        try:
            count = 0
            while True:
                count += 1
                #Calculate the centroid of the new cluster
                centsnew = []
                for one in clustlist:
                    if len(one) > 0:
                        m = [1.0]*len(one)
                        #cop = numpy.dot(m, one) / sum(m)
                        cop = [sum([es for es in one]) / sum(m)]
                        centsnew.append(cop)
                    else:
                        centsnew.append(one)
                #Assign individuals to new centroids
                clustlistn=[[] for i in range(Optimizer.nindiv)]
                for one in pts:
                    ds=[]
                    for i in range(len(centsnew)):
                        d = abs(one - cents[i])
                        ds.append([d,i])
                    ds.sort()
                    loc = ds[0]
                    clustlistn[loc[1]].append(one)
                if clustlistn==clustlist:
                    clustlist = clustlistn
                    centroids = centsnew
                    attemptcount=-10
                    passflag=True
                    break
                else:
                    clustlist = clustlistn
        except:
            attemptcount-=1
            STR+='WARNING: PREDATOR: K-means cluster difficulty\n'
            print 'WARNING: PREDATOR: K-means cluster difficulty'
    if passflag==True:
        #Select one individual from each cluster
        news = []
        for one in clustlist:
            #Random
            news.append(random.choice(one))
            #Minimum energy
            #one.sort()
            #news.append(one[0])

        indices = []
        for one in news:
            for i in range(len(pts)):
                if pts[i]==one:
                    indices.append(i)
                    break
        npop = [pop[ind] for ind in indices]
    else:
        STR+='WARNING: PREDATOR: Population converging. Unable to find nindiv distinct clusters\n'
        npop,str1 = Mutation_Dups(pop,Optimizer)
        STR+=str1
        #npop=pop
        STR+='WARNING: PREDATOR: Chose new population based on Mutation-Dups predator\n'
    pop = get_best(npop,len(npop))
    return pop