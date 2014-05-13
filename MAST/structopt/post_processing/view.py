import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy

def view(cluster, colorlist=[], nints=0):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	syms = list(set(cluster.get_chemical_symbols()))
	size=1
	for sym in syms:
		if sym in [s for s,c in colorlist]:
			colors = colorlist[[i for i,c in enumerate(colorlist) if c[0]==sym][0]][1]
		else:
			colors = numpy.random.rand(3,1)
			colorlist.append([sym,colors])
		xs = []
		ys = []
		zs = []
		for atm in cluster:
			if atm.index >= nints:
				if atm.symbol == sym:
					xs.append(atm.position[0])
					ys.append(atm.position[1])
					zs.append(atm.position[2])
					size = atm.mass**2
		ax.scatter(xs,ys,zs,s=size,c=colors)
	if 'defect' in [s for s,c in colorlist]:
		colors = colorlist[[i for i,c in enumerate(colorlist) if c[0]=='defect'][0]][1]
	#else:
	#	colors = numpy.random.rand(3,1)
	#	colorlist.append(['defect',colors])
	for sym in syms:
		xs = []
		ys = []
		zs = []
		for atm in cluster:
			if atm.index < nints:
				if atm.symbol == sym:
					xs.append(atm.position[0])
					ys.append(atm.position[1])
					zs.append(atm.position[2])
					size = atm.mass**2
		ax.scatter(xs,ys,zs,s=size,c=colors)
	ax.set_xlabel('X Position')
	ax.set_ylabel('Y Position')
	ax.set_zlabel('Z Position')
	plt.show()
	return colorlist