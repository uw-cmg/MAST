from matplotlib import pyplot as plt
from matplotlib import ticker
from MAST.structopt.tools.StemCalc import ConvStem
from MAST.structopt import Optimizer
from MAST.structopt import inp_out
import  numpy as np

aber=[[0,0],[0,0],[22.56,-20.1],[22.08,-7.5],[0.1198,0],[0.9018,-170.1],[0.04964,20.9],[28.43,-120.6],[11.84,153.8],[8.456,76.1],[0.622,0],[2.811,-125.5]]
autostemparameters={'Electron energy': 200,'Spherical aberration': 1.4,
'Defocus':0,'Aperture semiangle': 24.5,'Source size': 0.882,
'Slice size': 25.0,'Pixels':976,'Chromatic aberration Coefficient':1.4,'Delta E':0.73,
'aber':aber,'Scale Factor':0.00570113}
autostemparameters['Grid_sim2exp'] = 1
A = ConvStem(parameters=autostemparameters,calc_exp=False)

nk = autostemparameters['Pixels']
A.psf = np.empty([nk,nk],dtype=float)
fileobj = open('PSF.txt', 'r')
lines = fileobj.readlines()
for x in range(0,nk):
   A.psf[x] = lines[x].split()
fileobj.close()
#plt.imshow(A.psf)
#plt.show()

Au=inp_out.read_xyz('Output-rank0/indiv00.xyz',-1)
#Au=inp_out.read_xyz('STEM_ref',0)
imAu = A.get_image(A.psf, Au, autostemparameters['Slice size'], autostemparameters['Pixels'])
IPlot = imAu.T
plt.imshow(IPlot,origin="lower",cmap = plt.get_cmap('hot'))
cb=plt.colorbar()
tick_locator = ticker.MaxNLocator(nbins=7)
cb.locator = tick_locator
cb.update_ticks()
for t in cb.ax.get_yticklabels():
     t.set_fontsize(24)
plt.axis('off')
plt.show()

