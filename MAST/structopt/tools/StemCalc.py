import numpy
import os
import sys
import copy
import random
from tempfile import mkdtemp
import math
import cmath
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm  as cm
from PIL import Image, ImageChops
from ase.io import write
from ase import Atom, Atoms
import scipy.interpolate
import scipy.special
import scipy.misc
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
from scipy import optimize

class ConvStem():
    def __init__(self, parameters={}, mode='Simple', tmp_dir=None, keep_files=True, 
        calc_exp=True, debug=False):
        """Inputs:
        parameters - dictionary of parameters for ConvStem
        temp-dir - directory string to store ConvStem Files
        keep-Files - True/False for whether or not to keep files after run"""
        #Initialize parameters
        self.parameters=parameters
        self.keep_files=keep_files
        self.mode=mode
        self.debug=debug
        #Initialize directory for storing data
        if tmp_dir is None:
            self.tmp_dir = mkdtemp(prefix='ConvStem-')
        else:
            self.tmp_dir=os.path.realpath(tmp_dir)
            if not os.path.isdir(self.tmp_dir):
                os.mkdir(self.tmp_dir, 0755)
        #Get the probe function based on input parameters
        self.psf = get_probe_function(self.parameters)
        if calc_exp:
            if self.mode=='Simple':
                self.expfun=self.calculate_simp_function(self.parameters['Exp_Image'])
            elif self.mode=='Complex':
                self.expfun=self.calculate_comp_function(self.parameters['Exp_Image'])
            else:
                raise RuntimeError('Error: must specify either Simple or Complex mode')
            if self.keep_files:
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.imshow(self.expfun,cmap=cm.hot)
                ax.set_frame_on(False)
                ax.set_xticks([])
                ax.set_yticks([])
                plt.axis('off')
                ax.set_aspect('equal')
                figname = 'Experimental_Image.png'
                fig.savefig(self.tmp_dir+'/'+figname, transparent=True, bbox_inches='tight',pad_inches=0)
            
    def calculate_simp_function(self,image):
        #Load image file and convert image to grayscale
        #ima=Image.open(image).convert('LA')
        #a=numpy.array(ima)
        #a=scipy.misc.imread(image)
        #Get data array of image values
        #data=numpy.zeros((len(a),len(a[0])))
        #for i in range(len(a)):
        #    for j in range(len(a[0])):
        #        data[i][j]=a[i][j][0]
        #return data
        return image
    
    def calculate_comp_function(self,imagefile):
        #Function to calculate the overall gaussian function of an ConvStem simulation
        
        #Load image file and convert image to grayscale
        img=Image.open(imagefile).convert('LA')
        
        #Add ImageChops to account for strange pbc issue
        #imoff=ImageChops.offset(img,img.size[0]/2,img.size[1]/2)
        #a=numpy.array(imoff)
        
        #Get data array of image values
        data=numpy.zeros((len(a),len(a[0])))
        for i in range(len(a)):
            for j in range(len(a[0])):
                data[i][j]=a[i][j]
        
        #Calculate the atom column positions
        positions=get_atom_pos(self,data)
        x=positions[0]
        y=positions[1]
        
        #Apply rotation/translation
        mvs=[]
        for i in range(len(x)):
            mvs.append(data[y[i],x[i]])
        top2=sorted(enumerate(mvs), key=lambda a:a[1], reverse=True)[:2]
        #Set new origin at max value
        #Translation vector becomes:
        trans=[x[top2[0][0]],y[top2[0][0]]]
        #Find rotation angle alpha
        slope=float(y[top2[0][0]]-y[top2[1][0]])/float(x[top2[0][0]]-x[top2[1][0]])
        alpha=math.atan(abs(slope))
        #Adjust for quadrant
        if x[top2[1][0]] < x[top2[0][0]]:
            if y[top2[1][0]] > y[top2[0][0]]:
                #quadrant 2
                alpha=alpha+math.pi/2.0
            else:
                #quadrant 3
                alpha=alpha+math.pi
        else:
            if y[top2[1][0]] < y[top2[0][0]]:
                #quadrant 4
                alpha=alpha+math.pi*3.0/2.0
        xs=numpy.arange(len(data))
        xys=[]
        z=[]
        for i in range(len(xs)):
            for j in range(len(xs)):
                ara=numpy.array([xs[i], xs[j], 1])
                xys.append(ara)
                z.append(data[j][i])
        transfmat=numpy.array([[math.cos(alpha), -1.0*math.sin(alpha), -1.0*x[top2[0][0]]],[math.sin(alpha),math.cos(alpha),-1.0*y[top2[0][0]]],[0.0,0.0,1.0]])
        nx=[]
        ny=[]
        for one in xys:
            ara=numpy.dot(transfmat,one)
            nx.append(ara[0])
            ny.append(ara[1])
        
        #Calculate a gaussian intensity function for each atom
        #gauss=[]
        #for one in positions:
        #    gfun=get_gaussian(one)
        #    gauss.append(gfun)
        #Merge gaussian functions
        #gfuntot=merge_gauss(gauss)
        
        gfuntot=scipy.interpolate.RectBivariateSpline(x,y,z)
        
        return gfuntot
    
    def compare_functions(self, expfun, simfun):
        """Function compares two matrices and calculates chisq.  Matrices must be same size"""
        chi=[]
        for i in range(len(expfun)):
            for j in range(len(expfun[0])):
                c=(simfun[i][j]-expfun[i][j])**2/expfun[i][j]
                chi.append(c)
        chisq=sum(chi)
        return chisq
    
    def run(self,atms):
        #Get the image from ConvStem
        simfun=self.get_image(self.psf,atms,self.parameters['Slice size'],self.parameters['Pixels'])
        #Calculate the function for the image
        #if self.mode=='Simple':
        #    simfun=self.calculate_simp_function(image)
        #elif self.mode=='Complex':
        #    simfun=self.calculate_comp_function(image)
        chisq=self.compare_functions(self.expfun,simfun)
        if self.keep_files:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            #Uncomment the following lines to save image as borderless png
            ax.imshow(simfun,cmap=cm.hot)
            ax.set_frame_on(False)
            ax.set_xticks([])
            ax.set_yticks([])
            plt.axis('off')
            ax.set_aspect('equal')
            fignamep = 'SimImage-'
            #Label image with next available number
            x=0
            while True:
                figname = fignamep+repr(x)+'.png'
                if not os.path.exists(self.tmp_dir+'/'+figname):
                    break
                else:
                    x+=1
            fig.savefig(self.tmp_dir+'/'+figname, transparent=True, bbox_inches='tight',pad_inches=0)
        return chisq
    
    def get_image(self, psf, atms, rmax, nx, scalefactor=1.0):
        """Function to get image based point spread function and atoms
        rmax=Size of slice in Angstoms
        nx=number of pixels in slice"""
        #dr =  rmax / nx
        #ry = numpy.arange(0.0,rmax,dr)
        pot = stempot(rmax,rmax,len(psf),len(psf[0]),atms,scalefactor)
        #pot = stempot(xyz, rmax, rmax, nk, nk)
        potft = numpy.fft.fft2(pot)
        psf2d = numpy.fft.fft2(psf)
        potm = potft*psf2d
        zcon_im = numpy.fft.ifft2(potm)
        zcon_r = numpy.real(zcon_im)
        #zcon_s = numpy.fft.fftshift(zcon_r)
        zcon_s = zcon_r
        #grid = zcon_s.reshape((len(ry), len(ry)))
        #SetScale/P x DimOffSet(zcon_im, 0), DimDelta(zcon_im, 0), "", zcon_s
        #SetScale/P y DimOffSet(zcon_im, 1), DimDelta(zcon_im, 1), "", zcon_s
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)
    #     #Uncomment the following lines to save image as borderless png
    #     ax.imshow(grid,cmap=cm.hot)
    #     ax.set_frame_on(False)
    #     ax.set_xticks([])
    #     ax.set_yticks([])
    #     plt.axis('off')
    #     ax.set_aspect('equal')
    #     fignamep = 'SimImage-'
    #     #Label image with next available number
    #     x=0
    #     while True:
    #         figname = fignamep+repr(x)+'.png'
    #         if not os.path.exists(self.tmp_dir+'/'+figname):
    #             break
    #         else:
    #             x+=1
    #     fig.savefig(self.tmp_dir+'/'+figname, transparent=True, bbox_inches='tight',pad_inches=0)
        return zcon_s
    
    def get_atom_pos(self, data):
		"""Function to identify the location of the atom columns
		Inputs are:
			data = 2D matrix with imagefile value
			self.parameters:
				neighborhood_size = size of square stamp to search array
				threshold = difference for change"""
		
		# #Get data array of image values
		# data=list(ima.getdata())
		# imdat=[value for (value,max) in data]
		# width, height = ima.size
		# holder=numpy.zeros((height,width))
		# for h in height:
			# for w in width:
				# holder[h][w]=data[w+h*width]
		
		#Set values for peak guides
		if 'neighborhood_size' in self.parameters:
			neighborhood_size = self.parameters['neighborhood_size']
		else:
			neighborhood_size = 30
		if 'threshold' in self.parameters:
			threshold = self.parameters['threshold']
		else:
			threshold = 30
		
		#Use filters to calculate peaks
		data_max = filters.maximum_filter(data, neighborhood_size)
		maxima = (data == data_max)
		data_min = filters.minimum_filter(data, neighborhood_size)
		diff = ((data_max - data_min) > threshold)
		maxima[diff == 0] = 0

		labeled, num_objects = ndimage.label(maxima)
		slices = ndimage.find_objects(labeled)
		x, y = [], []
		for dy,dx in slices:
			x_center = (dx.start + dx.stop - 1)/2
			x.append(x_center)
			y_center = (dy.start + dy.stop - 1)/2    
			y.append(y_center)
		
		#Plot peaks points overlayed on image
		#plt.imshow(data)
		#plt.scatter(x,y)
		#plt.show()
		
		#positions=[[x[i],y[i]] for i in range(len(x))]
		posiitons=[x,y]
		
		return positions
    
def get_probe_function(parameters):
    """Function to get the probe function based on input parameters
    kev=Electron energy in keV
    ap = Objective aperture semiangle in mrad
    Cc = Chromatic aberration coefficient in mm
    dE = Delta E in eV
    ds=Source size in Angstroms
    rmax=Size of slice in Angstroms
    Cs=Spherical aberration Cs in mm
    df=Defocus in Angstroms
    nx=number of pixels in slice"""
    
    keV = parameters['Electron energy']
    ap = parameters['Aperture semiangle']
    try:
        Cc = parameters['Chromatic aberration Coefficient']
    except:
        Cc = 0.0
    try:
        dE = parameters['Delta E']
    except:
        dE = 0.0
    #Cs=parameters['Spherical aberration']
    #df=parameters['Defocus']
    ds=parameters['Source size']
    rmax=parameters['Slice size']
    #nx=parameters['Pixels']
    try:
        aber=parameters['aber']
    except:
        aber=[[0,0] for i in range(12)]
    
    nk = math.floor(40*(0.001*ap / wavlen(keV))*rmax)
    if math.fmod(nk, 2):
        nk += 1
    
    #Get probe function
    if Cc == 0.0 and dE == 0.0:
        psf2D = STEMPSF2DCoh(aber, keV, ap, nk)
        #Redimension/D/C probe2DCoh #Converts wave to double precision complex wave
        #wave/c psf2D = $"probe2DCoh"
    else:
        psf2D = STEMPSF2DIncoh(aber, keV, Cc, dE, ds, ap, nk)
        #Redimension/D/C probe2DIncoh #Converts wave to double precision complex wave
        #wave/c psf2D = $"probe2DIncoh"
    #Shift minimum to zero
    mn = numpy.minimum.reduce(numpy.minimum.reduce(psf2D))
    for i in range(len(psf2D)):
        for j in range(len(psf2D[0])):
            psf2D[i][j] += abs(mn)
    
    return psf2D

def STEMPSF2DCoh(aber, keV, ap, nk):
    kmax = 0.001*ap / wavlen(keV)
    phase = ChiPhase2D(aber, keV, ap, nk/10)
    xp=numpy.linspace(-kmax*2.0,(kmax*2.0),nk/10)
    yp=numpy.linspace(-kmax*2.0,(kmax*2.0),nk/10)
    phasef = scipy.interpolate.RectBivariateSpline(xp,yp,phase)
    
    # need to pad the phase with zeros here.
    probe2DCoh = numpy.zeros((nk,nk),dtype=complex)
    x=numpy.linspace(-20*kmax,40*kmax,nk)
    y=numpy.linspace(-20*kmax,40*kmax,nk)
    
    #SetScale/I x -20*kmax, 40*kmax, "", probe2DCoh
    #SetScale/I y -20*kmax, 40*kmax, "", probe2DCoh
    for i in range(len(x)):
        for j in range(len(y)):
            if x[i]**2 + y[j]**2 < kmax**2:
                probe2DCoh[i][j] = complex(1, phasef(x[i],y[j]))
            else:
                probe2DCoh[i][j] = complex(0, 0)
    #probe2DCoh = ( (x^2 + y^2 < kmax^2) ?  p2rect(cmplx(1, phase(x)(y)))  : cmplx(0, 0))

    p2dcf = numpy.fft.fft2(probe2DCoh)
    for i in range(len(x)):
        for j in range(len(y)):
            probe2DCoh[i][j] = complex((abs(p2dcf[i][j]))**2, 0)
    probe2DCoh = numpy.real(probe2DCoh)    
    probe2DCoh /= sum(sum(probe2DCoh))
    
    return probe2DCoh

def STEMPSF2DIncoh(aber, keV, Cc, dE, ds, ap, nk):
    kmax = 0.001*ap / wavlen(keV)
    defocus_distribution = ChromaticDefocusDistribution(Cc, dE, keV, ap)
    start_df = aber[0][0]
    df_range = 2.5*(Cc*dE /(1e3*keV)*( (1+keV/511)/(1+keV/1022) ))
    x = numpy.linspace(-df_range,df_range, num=len(defocus_distribution))
    for i in range(len(defocus_distribution)):
        aber[0][0] = start_df + x[i]
        probe2DCoh = STEMPSF2DCoh(aber, keV, ap, nk)
        if i==0:
            probe2DIncoh = copy.deepcopy(probe2DCoh)
            probe2DIncoh *= defocus_distribution[0]
        else:
            for j in range(len(probe2DIncoh)):
                for k in range(len(probe2DIncoh[0])):
                    probe2DIncoh[j][k] += defocus_distribution[i]*probe2DCoh[j][k]
    
    probe2DIncoh /= sum(defocus_distribution)

    aber[0][0] = start_df

    if ds != 0:
        fds = ds / (2*(2*math.log(2))**0.5)    # real-space standard deviation for Gaussian with FWHM ds
        fds = 1/(2*math.pi*fds)    # FT standard deviation
        #Redimension/C probe2DIncoh
        #wave/C probe2DSS = $"probe2DIncoh"
        probe2DSS = numpy.fft.fft2(probe2DIncoh)
        g = Gauss(0.0, 0.0, fds, fds)#, len(probe2DSS), len(probe2DSS[0]))
        xs = numpy.linspace(-1.0, 1.0, num=len(probe2DSS))
        for j in range(len(probe2DSS)):
                for k in range(len(probe2DSS[0])):
                    probe2DSS[j][k]*=complex(g(xs[j],xs[k]),0.0)
        #probe2DSS *= complex(Gauss(x, 0.0, fds, y, 0.0, fds), 0.0)
        probe2DSS = numpy.fft.ifft2(probe2DSS)
        probe2DSS = numpy.real(probe2DSS)
        probe2DSS = numpy.fft.fftshift(probe2DSS)
        #probe2DSS = cmplx(sqrt(magsqr(probe2DSS)), 0)
        probe2DSS /= sum(sum(probe2DSS))
        probe2DIncoh = copy.deepcopy(probe2DSS)
        #SetScale/p x dimoffset(probe2dcoh, 0), dimdelta(probe2dcoh, 0), "", probe2dIncoh
        #SetScale/P y dimoffset(probe2dcoh, 1), dimdelta(probe2dcoh, 1), "", probe2dIncoh
    
    return probe2DIncoh

def ChiPhase2D(aber, keV, ap, nk):
    aber = SwitchToAngstroms(aber)
    wl = wavlen(keV)
    kmax = (0.001*ap / wl)  # maximum k through aperture
    x=numpy.linspace(-kmax*2.0,(kmax*2.0),nk)
    y=numpy.linspace(-kmax*2.0,(kmax*2.0),nk)
    astack = [numpy.zeros((len(x),len(y))) for one in range(12)]
    
    # Evaluate the phase shifts from the various aberrations
    for i in range(len(x)):
        for j in range(len(y)):
            #C1, defocus
            astack[0][i][j] = (1.0/2.0)*wl*aber[0][0]*(x[i]**2.0 + y[j]**2.0)
            #A1, 2-fold astigmatism
            astack[1][i][j] = (1.0/2.0)*wl*aber[1][0]*(x[i]**2.0 - y[j]**2.0)
            #A2, 3-fold astigmatism
            astack[2][i][j] = (1.0/3.0)*wl**2.0*aber[2][0]*(x[i]**3.0 - 3.0*x[i]*y[j]**2.0)
            #B2, axial coma
            astack[3][i][j] = wl**2.0*aber[3][0]*(x[i]**3.0 + x[i]*y[j]**2.0)
            #C3, primary spherical aberration
            astack[4][i][j] = (1.0/4.0)*wl**3.0*aber[4][0]*(x[i]**4.0 + 2.0*x[i]**2.0*y[j]**2.0 + y[j]**4.0)
            #A3, 4-fold astigmatism
            astack[5][i][j] = (1.0/4.0)*wl**3.0*aber[5][0]*(x[i]**4.0 - 6.0*x[i]**2.0*y[j]**2.0 + y[j]**4.0)
            #S3, star aberration
            astack[6][i][j] = wl**3.0*aber[6][0]*(x[i]**4.0 - y[j]**4.0)
            #A4, 5-fold astigmatism
            astack[7][i][j] = (1.0/5.0)*wl**4.0*aber[7][0]*(x[i]**5.0 - 10.0*x[i]**3.0*y[j]**2.0 + 5.0*x[i]*y[j]**4.0)
            #D4, 3-lobe aberration
            astack[8][i][j] = wl**4.0*aber[8][0]*(x[i]**5.0 - 2.0*x[i]**3.0*y[j]**2.0 - 3.0*x[i]*y[j]**4.0)
            #B4, axial coma
            astack[9][i][j] = wl**4.0*aber[9][0]*(x[i]**5.0 + 2.0*x[i]**3.0*y[j]**2.0 + x[i]*y[j]**4.0)
            #C5, 5th order spherical aberration
            astack[10][i][j] = (1.0/6.0)*wl**5.0*aber[10][0]*(x[i]**6.0 + 3.0*x[i]**4.0*y[j]**2.0 + 3.0*x[i]**2.0*y[j]**4.0 + y[j]**6.0)
            #A5, 5th order spherical aberration
            astack[11][i][j] = (1.0/6.0)*wl**5.0*aber[11][0]*(x[i]**6.0 - 15.0*x[i]**4.0*y[j]**2.0 + 15.0*x[i]**2.0*y[j]**4.0 - y[j]**6.0)
    
    #Set minimum to zero
    #nnastack = [numpy.zeros((len(x),len(y))) for one in range(12)] 
    # for i in range(len(x)):
#         for j in range(len(y)):
#             for k in range(len(astack)):
#                 if astack[k][i][j] < 0.0:
#                     astack[k][i][j]=0.0
    #Shift minimum to zero
    for k in range(len(astack)):
        mn = numpy.minimum.reduce(numpy.minimum.reduce(astack[k]))
        for i in range(len(x)):
            for j in range(len(y)):
                astack[k][i][j]+= abs(mn)
    
    # rotate the phase shifts of the non-centrosymmetric aberrations
    for i in range(12):
        if aber[i][1] != 0.0:
            astack[i] = scipy.ndimage.rotate(astack[i],-aber[i][1])
            #ImageTransform/P=(i) getplane astack
            #wave aphase = $"M_ImagePlane"
            #ImageRotate/E=0/O/A=(aber[i][1]) aphase
            #SetScale/P x -DimSize(aphase, 0)*DimDelta(astack, 0) / 2, DimDelta(astack, 0), "", aphase
            #SetScale/P y -DimSize(aphase, 1)*DimDelta(astack, 1) / 2, DimDelta(astack, 0), "", aphase
            #astack[][][i] = aphase(x)(y)

    # sum all the aberration contributions
    fsum = numpy.zeros((len(x),len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            fsum[i][j] = sum([one[i][j] for one in astack])*2*math.pi
    #MatrixOp/O phase = 2*Pi*sumbeams(astack)
    #SetScale/I x -2*kmax, 2*kmax, "", phase
    #SetScale/I y -2*kmax, 2*kmax, "", phase
    
    return fsum

def ChromaticDefocusDistribution(Cc, dE, keV, ap):
    Cc *= 1e7  # mm to Angstroms
    df_phase_max = 2*math.pi / 50  # maximum phase step at the aperture edge due to chromatic aberration
    kmax = 0.001*ap/wavlen(keV)
    # defocus range and form from Reimer
    H = (Cc*dE /(1e3*keV)*( (1+keV/511)/(1+keV/1022) ))
    N = (2*(math.log(2)**0.5) / ((math.pi**0.5)*H))
    df_range = 2.5*H
    ndf = math.ceil(df_range * wavlen(keV) * kmax**2 / df_phase_max)
    if ndf < 31:
        ndf = 31
    
    #ndf = (ndf < 31 ? 31 : ndf)
    if not math.fmod(ndf,2):
        ndf+=1
    
    #ndf = (!mod(ndf, 2) ? ndf+1 : ndf)
    defocus_distribution = numpy.linspace(-df_range, df_range, num = ndf)
    #Make/O/N=(ndf) defocus_distribution
    #SetScale/I x -df_range, df_range, "", defocus_distribution
    for x in range(len(defocus_distribution)):
        defocus_distribution[x] = N*math.exp( -math.log(2) * (2*defocus_distribution[x] / H)**2 )
    
    return defocus_distribution

def stempot(xmax,ymax,nx,ny,atms,scalefactor):
    """function V = stempot(xmax,ymax,nx,ny,potfile)
    %  STEMPOT Generate a Projected Potential
    %  inputs xmax, ymax are the size of the slice in angstroms
    %  nx,ny are the number of pixels in the x and y directions """
    #zed=2 for rutherford scattering of the nucleus, less for screening
    zed = 1.7

    ix = numpy.arange(1.0,nx)
    iy = numpy.arange(1.0,ny)
    dx = xmax/nx
    dy = ymax/ny
    rx = numpy.arange(0,xmax-dx,dx)
    ry = numpy.arange(0,ymax-dy,dy)

    Zatom = atms.get_atomic_numbers()
    #translate atoms such that the center of mass is in the center of the computational cell
    com = atms.get_center_of_mass()
    cop = xmax/2.0
    trans = [cop-i for i in com]
    atms.translate(trans)
    positions=atms.get_positions()
    ax=[]
    ay=[]
    az=[]
    for o,t,h in positions:
        ax.append(o)
        ay.append(t)
        az.append(h)
    ax = numpy.array(ax)
    ay = numpy.array(ay)
    az = numpy.array(az)
    amax = len(Zatom)

    #find boundaries of slice
    axmin = min(ax)
    axmax = max(ax)
    aymin = min(ay)
    aymax = max(ay)
    
    V= numpy.zeros((nx,ny))

    #map x and y coords of the atoms to the nearest grid points
    #A fraction of the atom must be assigned to the closest gridpoints
    #to avoid sum and difference frequencies appearing in the image
    #grid point to the left of the atom
    iax = numpy.array([math.floor(axi/dx)+1 for axi in ax])
    #create periodic boundary conditions
    ibx=numpy.array([math.fmod(iaxi,nx)+1 for iaxi in iax])
    #fraction of atom at iax
    fax = numpy.array([1-math.fmod((axi/dx),1 ) for axi in ax]) 
    #grid point above the atom
    iay = numpy.array([math.floor(ayi/dy)+1 for ayi in ay])
    #create periodic boundary conditions
    iby = numpy.array([math.fmod(iayi,ny)+1 for iayi in iay])
    #fraction of atom at iay 
    fay = numpy.array([1-math.fmod((ayi/dy),1 ) for ayi in ay])
    
    #Add each atom to the potential grid
    V1 = numpy.array([fax[i] * fay[i] * (Zatom[i]**zed) for i in range(len(fax))])
    V2 = numpy.array([(1-fax[i]) * fay[i] * (Zatom[i]**zed) for i in range(len(fax))])
    V3 = numpy.array([fax[i] * (1-fay[i]) * (Zatom[i]**zed) for i in range(len(fax))])
    V4 = numpy.array([(1-fax[i]) * (1-fay[i]) * (Zatom[i]**zed) for i in range(len(fax))])
#     V1 = numpy.array([fax[i] * fay[i] * scalefactor for i in range(len(fax))])
#     V2 = numpy.array([(1-fax[i]) * fay[i] * scalefactor for i in range(len(fax))])
#     V3 = numpy.array([fax[i] * (1-fay[i]) * scalefactor for i in range(len(fax))])
#     V4 = numpy.array([(1-fax[i]) * (1-fay[i]) * scalefactor for i in range(len(fax))])

    for j in range(amax):
       V[iax[j],iay[j]] += V1[j]
       V[ibx[j],iay[j]] += V2[j]
       V[iax[j],iby[j]] += V3[j]
       V[ibx[j],iby[j]] += V4[j]
    
    return V

def wavlen(keV):
    #calculate electron wavelength in Angstroms given the energy in keV.
    wav = 12.3986 / ( (2*511.0 + keV) * keV)**0.5
    return wav

def SwitchToAngstroms(aber):
    '''Function to switch aberrations wave into Angstroms instead of natural units'''
    naber = [[0.0,0.0] for i in range(12)]
    #C1, A1, A3, B2, all start in nm
    for i in range(0,4):
        naber[i][0] = 10*aber[i][0]
    #C3, A3, S3, A4, D4, B4 start in um
    for i in range(4,10):
        naber[i][0] = 1e4*aber[i][0]
    #C5, A5 in mm
    for i in range(10,12):
        naber[i][0] = 1e7*aber[i][0]
    
    return naber

def Gauss(center_x, center_y, width_x, width_y, height=1.0):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*math.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)

def moments(data):
	"""Returns (height, x, y, width_x, width_y)
		the gaussian parameters of a 2D distribution by calculating its moments """
	total = data.sum()
	X, Y = indices(data.shape)
	x = (X*data).sum()/total
	y = (Y*data).sum()/total
	col = data[:, int(y)]
	width_x = sqrt(abs((arange(col.size)-y)**2*col).sum()/col.sum())
	row = data[int(x), :]
	width_y = sqrt(abs((arange(row.size)-x)**2*row).sum()/row.sum())
	height = data.max()
	return x, y, width_x, width_y, height

def fitgaussian(data):
	"""Returns (height, x, y, width_x, width_y)
		the gaussian parameters of a 2D distribution found by a fit"""
	params = moments(data)
	errorfunction = lambda p: ravel(Gauss(*p)(*indices(data.shape)) - data)
	p, success = optimize.leastsq(errorfunction, params)
	return p

def find_stem_coeff(Optimizer, indiv):
    from MAST.structopt.tools.eval_energy import eval_energy
    outs = eval_energy([Optimizer,indiv])
    indiv.energy = outs[0]
    stro=outs[3]
    if Optimizer.structure == 'Defect' or Optimizer.structure=='Surface':
        indiv.bulki = outs[1]
    chisq = Optimizer.stemcalc.run(indiv[0])
    aflag=True
    alpha = 1.0
    while True:
        value=alpha*chisq
        div=abs(indiv.energy)/value
        if div <1:
            alpha*=0.1
        elif div >10:
            alpha*=10
        else:
            break
    indiv.fitness=indiv.energy+alpha*chisq
    return alpha, indiv
