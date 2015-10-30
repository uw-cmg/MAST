import matplotlib.pyplot as plt

def read_genealogy(filename):
    """Function to read genealogy and generate individual tree
    NEEDS DEVELOPMENT
    """
    f = open(filename,'r')
    line0=f.readline()
    indhist=[[value] for value in line0.split()]
    for line in f.readlines():
        i=0
        for value in line.split():
            indhist[i].append(value)
            i+=1
    f.close()
    numinds=len(indhist)
    gens=len(indhist[0])
    xs=[1000/numinds*(n+1) for n in range(numinds)]
    ys=[1000-1000.0/float(gens)*(n+1) for n in range(gens)]
    fig = plt.figure()
    #figsize=10 #max(numinds,gens)
    fig.set_size_inches(numinds, gens, forward=True)
    ax1=fig.add_subplot(111)
    ax1.set_frame_on(False)
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.axis('off')
    for i in range(numinds):
        prev=None
        for j in range(gens):
            if indhist[i][j] !=prev:
                ax1.text(xs[i],ys[j],indhist[i][j])
                prev=indhist[i][j]
    #plt.xlim([0,1000])
    #plt.ylim([0,1000])
    lines=[]
    for i in range(numinds):
        for j in range(gens):
            if ')' in indhist[i][j]:
                a=indhist[i][j].split('(')[1].split(')')
                ind=[int(one) for one in a[0].split('+')]
                if 'm' in indhist[i][j]:
                    lines.append([xs[ind[0]],ys[j-1],xs[i],ys[j],'b'])
                    lines.append([xs[ind[1]],ys[j-1],xs[i],ys[j],'b'])
                else:
                    lines.append([xs[ind[0]],ys[j-1],xs[i],ys[j],'k'])
                    lines.append([xs[ind[1]],ys[j-1],xs[i],ys[j],'k'])
            elif 'm' in indhist[i][j]:
                #TTM 20151030 retab on the document to change tabs to '    '.
                a=indhist[i][j].split('m')
                ind=int(float(a[0])) # TTM 20151030 float can convert a string, so make it a float first.
                lines.append([xs[ind],ys[j-1],xs[i],ys[j],'b'])

    for one in lines:
        ax1.plot([one[0], one[2]], [one[1], one[3]], color=one[4])

    plt.savefig('geneologytree-'+filename+'.png')

