import matplotlib.pyplot as plt
import os
import numpy
from MAST.structopt.post_processing.read_genealogy import read_genealogy
import logging

def read_output(folder,genealogytree=False,natoms=None, loggername=None):
    """Function to analyze structural optimization output.
    Inputs:
        folder = folder containing output from structural optimization
        genealogytree = Boolean for whether or not to construct genealogy tree from data
            Default = False
        natoms = number of atoms in an individual 
            Default = False
        loggername = name of log file to write to
    Outputs:
        Plot-All_Energies-XXX.png
        Plot-Peratom-All_Energies-XXX.png
        Plot-Summary-XXX.png
        Plot-time-Summary-XXX.png
        Plot-FPMin-FingerprintMin-XXX-genXX.png
        Plot-FpDist-Fingerprints-XXX-genXX.png
        Plot-SuccessfulMutations-XXX.png
        genealogytree-Genealogy-XXX.png
        
    """
    if loggername:
        #logger = initialize_logger(loggername)
        logger = logging.getLogger(Optimizer.loggername)
    files=os.listdir(folder)
    for filename in files:
        if 'All_Energies' in filename:
            if loggername:
                logger.info('Reading from All_Energies file.  Plotting all energies vs. generation')
            f=open(os.path.join(folder,filename),'r')
            line0=f.readline()
            indenergies=[[float(value)] for value in line0.split()]
            for line in f.readlines():
                i=0
                for value in line.split():
                    a=float(value)
                    if a != 0:
                        indenergies[i].append(float(value))
                    else:
                        indenergies[i].append(None)
                    i+=1
    
            f.close()
            generation=range(len(indenergies[0]))
            clist=[]
            for i in range(len(indenergies)):
                clist.append(numpy.random.rand(3,1))
    
            fig=plt.figure()
            ax1=fig.add_subplot(111)
            i=0
            for i in range(len(indenergies)):
                ax1.scatter(generation,indenergies[i],c=clist[i],label='Ind '+repr(i+1))
    
            handles, labels = ax1.get_legend_handles_labels()
            ax1.legend(handles, labels)
            box = ax1.get_position()
            ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            #ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.0),ncol=4, fancybox=True, shadow=True)
            plt.xlabel('Generation')
            plt.ylabel('Total Energy, eV')
            plt.title('Energy Evolution')
            plt.xlim([0,len(generation)])
            plt.savefig('Plot-{0}.png'.format(filename))
            indperatom=[]
            if natoms !=None:
                if loggername:
                    logger.info('Plotting energy per atom vs. generation')
                for ind in indenergies:
                    indperatom.append([en/natoms for en in ind])
        
                fig=plt.figure()
                ax1=fig.add_subplot(111)
                i=0
                for i in range(len(indperatom)):
                    ax1.scatter(generation,indperatom[i],c=clist[i],label='Ind '+repr(i+1))
        
                handles, labels = ax1.get_legend_handles_labels()
                ax1.legend(handles, labels)
                box = ax1.get_position()
                ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.xlabel('Generation')
                plt.ylabel('Energy/atom, eV/atom')
                plt.title('Energy Evolution')
                plt.xlim([0,len(generation)])
                plt.savefig('Plot-Peratom-{0}.png'.format(filename))
        if 'Summary' in filename:
            if 'StructureSummary' not in filename:
                if loggername:
                    logger.info('Plotting Fitness min, max and medium vs. generation')
                f=open(os.path.join(folder,filename),'r')
                line0=f.readline()
                line1=f.readline().split(':')
                bulkepa=float(line1[1])
                line2=f.readline().split()
                natoms=float(line2[1])
                labels=f.readline().split()
                data=dict((lab,[]) for lab in labels)
                for line in f.readlines():
                    i=0
                    for value in line.split(' ',6):
                        if i <6:
                            data[labels[i]].append(float(value))
                        else:
                            sp=value.split("'")
                            data[labels[i]].append(sp[1])
                        i+=1
    
                f.close()
                yl=[data['Fitmedium'][one]-data['Fitmin'][one] for one in range(len(data['Generation']))]
                yh=[data['Fitmax'][one]-data['Fitmedium'][one] for one in range(len(data['Generation']))]
                fig=plt.figure()
                ax1=fig.add_subplot(111)
                ax1.errorbar(data['Generation'],data['Fitmedium'],xerr=0,yerr=[yl,yh])
                plt.xlabel('Generation')
                plt.ylabel('Fitness')
                plt.xlim([1,len(data['Generation'])+1])
                plt.title('Fitness vs. Generation (Medium, Min, Max)')
                plt.ylim([min(data['Fitmin'])-1.0,max(data['Fitmin'])+10.0])
                plt.savefig('Plot-{0}.png'.format(filename))
                srtsl=data['time'][0].split()
                strtime=[float(one) for one in srtsl[3].split(':')]
                tclk=[0.0]
                for one in data['time']:
                    stl=one.split()
                    otime=[float(two) for two in stl[3].split(':')]
                    diff=[otime[i]-strtime[i] for i in range(len(otime))]
                    if diff[0]<0:
                        diff[0]+=24
                    for i in range(len(diff)):
                        if diff[i]<0:
                            diff[i-1]-=1.0
                            diff[i]+=60
                    tclk.append(diff[0]*3600+diff[1]*60+diff[2])
                    strtime=otime
                if loggername:
                    logger.info('Plotting computational time per generation')
                fig = plt.figure()
                ax1 = fig.add_subplot(111)
                ax1.plot(tclk)
                plt.xlabel('Generation')
                plt.ylabel('CPU Time (s)')
                plt.xlim([1,len(data['Generation'])+1])
                plt.title('Processing Time per Generation')
                totaltime=sum(tclk)/3600.0
                plt.text(0.5, 0.95,'Total Time for run = '+'%.3f' % round(totaltime, 3)+'hours',horizontalalignment='center',verticalalignment='top',transform = ax1.transAxes)
                plt.savefig('Plot-time-{0}.png'.format(filename))
        if 'FingerprintMin' in filename:
            if loggername:
                logger.info('Plotting minimum fingerprint functions by generation')
            f = open(os.path.join(folder,filename),'r')
            fpmin=[]
            enfp=[]
            i=0
            for line in f.readlines():
                if i==0:
                    bins=line.split(', ')
                    nbin=[]
                    for one in bins:
                        if '[' in one:
                            one=one[1::]
                        elif ']'in one:
                            ind=[i for i,some in list(enumerate(one)) if some==']']
                            one=one[0:ind[0]]
                        nbin.append(float(one))
                    fpmin.append(nbin)
                    i=1
                elif i==1:
                    enfp.append(float(line))
                    i=0
            f.close()
            for i in range(len(fpmin)):
                fig=plt.figure()
                ax=fig.add_subplot(111)
                ax.plot(fpmin[i])
                plt.xlabel('Distance, bins')
                plt.ylabel('Intensity')
                plt.title('Minimum Individual Fingerprint, Generation={0}'.format(i))
                plt.text(0.5, 0.95,'Structure Energy = '+repr(enfp[i]),horizontalalignment='center',verticalalignment='top',transform = ax.transAxes)
                if max(fpmin[i])>100:
                    plt.ylim([0,100])
                plt.savefig('Plot-FPMin-{0}-gen{1}.png'.format(filename,i))
        if 'Fingerprints' in filename:
            if loggername:
                logger.info('Plotting fingerprint distance vs. energy by generation')
            f = open(os.path.join(folder,filename),'r')
            fpds=[]
            ens=[]
            for line in f.readlines():
                fpg=[]
                eng=[]
                i=0
                for value in line.split():
                    if i==0:
                        fpg.append(float(value))
                        i=1
                    else:
                        eng.append(float(value))
                        i=0
                fpds.append(fpg)
                ens.append(eng)
    
            f.close()
    
            for i in range(len(fpds)):
                emin=min(ens[i])
                eng=[one-emin for one in ens[i]]
                fig=plt.figure()
                ax=fig.add_subplot(111)
                ax.scatter(fpds[i],eng)
                plt.xlabel('Fingerprint Cosine Distance from Minimum Energy Structure')
                plt.ylabel('Energy Difference, eV')
                plt.title('Fingerprint Distance vs. Energy, Generation={0}'.format(i))
                plt.savefig('Plot-FpDist-{0}-gen{1}.png'.format(filename,i))
        if 'Genealogy' in filename:
            if loggername:
                logger.info('Plotting bar plot of mutation success')
            f = open(os.path.join(folder,filename),'r')
            line0=f.readline()
            indhist=[[value] for value in line0.split()]
            for line in f.readlines():
                i=0
                for value in line.split():
                    indhist[i].append(value)
                    i+=1
            mutoptn = ['Lattice_Alteration_nn', 'Lattice_Alteration_rdrd', 'Lattice_Alteration_Group', 'Rotation_geo', 'ZP_Rotation', 'Random_Replacement']
            mutopt = ['LANN','LARD','LAG','GR','ZPR','RGR']
            count = [[0 for j in range(len(indhist[0]))] for i in mutopt]
            for i in range(len(indhist[0])):
                for j in range(len(indhist)):
                    for k in range(len(mutopt)):
                        if mutopt[k] in indhist[j][i]:
                            count[k][i]+=1
            generation=numpy.arange(len(indhist[0]))
            width=1.0/(len(count)+1)
            clist=[]
            for i in range(len(count)):
                clist.append(numpy.random.rand(3,1))

            fig=plt.figure()
            ax1=fig.add_subplot(111)
            for i in range(len(count)):
                ax1.bar(generation+width*i,count[i],width,color=clist[i],label=mutoptn[i])

            plt.xlabel('Generation')
            plt.ylabel('Successful Mutations')
            ax1.legend(loc='upper left')
            plt.title('Successful Mutations by Generation')
            plt.savefig('Plot-SuccessfulMutations-{0}.png'.format(filename))
            if genealogytree:
                if loggername:
                    logger.info('Plotting genealogy tree')
                read_genealogy(filename)
        if 'Bests-energies' in filename:
            if loggername:
                logger.info('Plotting Best-energies scatter')
            f=open(os.path.join(folder,filename),'r')
            benergies = []
            for line in f.readlines():
                benergies.append(float(line))    
            f.close()
            fig=plt.figure()
            ax1=fig.add_subplot(111)
            ax1.scatter(range(len(benergies)),benergies)
            plt.xlabel('Rank')
            plt.ylabel('Total Energy, eV')
            plt.title('Energy of Best Structures in Optimization')
            plt.savefig('Plot-{0}.png'.format(filename))
        else:
            pass

