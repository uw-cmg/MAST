#import matplotlib.pyplot as plt
#import matplotlib.cm  as cm
from MAST.structopt.tools.StemCalc import ConvStem
from MAST.structopt.generate.gen_pop_box import gen_pop_box as gpb
from MAST.structopt import Optimizer
import numpy

size = 10.0
Au1 = gpb([('Au',10,0,0),('Cu',10,0,0)],size)
aber=[[0,0],[0,0],[22.56,-20.1],[22.08,-7.5],[0.1198,0],[0.9018,-170.1],[0.04964,20.9],[28.43,-120.6],[11.84,153.8],[8.456,76.1],[0.622,0],[2.811,-125.5]]
autostemparameters={'Electron energy': 200,'Spherical aberration': 1.4,'Defocus':0,'Aperture semiangle': 24.5,'Source size': 0.882,'Slice size': 20.0,'Pixels':800,'aber':aber,'Scale Factor':0.00570113}
A = ConvStem(parameters=autostemparameters,calc_exp=False)
#fig = plt.figure()
#ax = fig.add_subplot(111)
#psf = numpy.fft.fftshift(A.psf)
#ax.imshow(psf,cmap=cm.hot)
#plt.show()
imAu1 = A.get_image(A.psf, Au1, autostemparameters['Slice size'], autostemparameters['Pixels'], 1)
#fig2 = plt.figure()
#ax2 = fig2.add_subplot(111)
#ax2.imshow(imAu1,cmap=cm.hot)
#plt.show()
autostemparameters['Exp_Image'] = imAu1
optparameters = {'structure': 'Cluster', 'optimizer_type': 'BH', 'size':10,'atomlist': [('Au', 10, 0, 0), ('Cu', 10, 0, 0)],'FIT_SCHEME':'STEM_Cost','STEM_Parameters':autostemparameters}
B = Optimizer(optparameters)
B.algorithm_serial()
