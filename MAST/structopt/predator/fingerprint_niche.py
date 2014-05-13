from MAST.structopt.fingerprinting import fingerprint_dist
from MAST.structopt.tools import get_best
from MAST.structopt.predator import mutation_dups
import random

def fingerprint_niche(pop, Optimizer):
    """Use a k-means clustering approach to identify individuals for new population"""
    STR = ''
    #Get coordinates for each individual
    pts = []
    for one in pop:
        fpd = fingerprint_dist(pop[0].fingerprint,one.fingerprint)
        pts.append([one.fitness - pop[0].fitness,fpd])

    for one in pts:
        if numpy.isnan(one[1]):
            one[1] = 0
    #Identify spatial bounds
    mins = min(pts)
    maxs = max(pts)

    attemptcount=10
    passflag=False
    while attemptcount > 0:
        #Select random initial centroid locations
        cents = []
        #for i in range(Optimizer.nindiv):
        #	cents.append([random.uniform(mins[0],maxs[0]),random.uniform(mins[1],maxs[1])])
        while len(cents) < Optimizer.nindiv:
            a = random.choice(pts)
            if a not in cents:
                cents.append(a)
        #Assign individual to closest centroid
        clustlist = [[] for i in range(Optimizer.nindiv)]
        for one in pts:
            ds = []
            for i in range(len(cents)):
                d = ((one[0] - cents[i][0])**2 + (one[1] - cents[i][1])**2)**0.5
                ds.append([d,i])
            ds = sorted(ds, key=lambda two: two[0])
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
                        cop = [sum([es for es,fps in one]) / sum(m),sum([fps for es,fps in one]) / sum(m)]
                        centsnew.append(cop)
                    else:
                        centsnew.append(one)
                #Assign individuals to new centroids
                clustlistn=[[] for i in range(Optimizer.nindiv)]
                for one in pts:
                    ds=[]
                    for i in range(len(centsnew)):
                        d = ((one[0] - centsnew[i][0])**2 + (one[1] - centsnew[i][1])**2)**0.5
                        ds.append([d,i])
                    ds = sorted(ds, key=lambda two: two[0])
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
            STR+= 'WARNING: PREDATOR: K-means cluster difficulty\n'
    if passflag==True:
        #Select one individual from each cluster
        news = []
        for one in clustlist:
            #Random
            #npop.append(random.choice(one))
            #Minimum energy
            opts = sorted(one, key=lambda two: two[0])
            news.append(opts[0])

        indices = []
        for one in news:
            for i in range(len(pts)):
                if pts[i]==one:
                    indices.append(i)
                    break
        npop = [pop[ind] for ind in indices]
    else:
        STR+='WARNING: PREDATOR: Population converging. Unable to find nindiv distinct clusters\n'
        npop = mutation_dups(pop,Optimizer)
        #npop = pop
        STR+='WARNING: PREDATOR: Chose new population based on Mutation-Dups predator\n'
    pop = get_best(npop,len(npop))
    return pop, STR