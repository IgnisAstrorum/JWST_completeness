#CODE TO PRODUCE A LIBRARY OF GALAXY TEMPLATES, THESE CAN THEN BE SCALED BY A MAGNITUDE AND ADDED TO AN IMAGE LATER

# coding: utf-8

# In[1]:

from multiprocessing import Pool
from joblib import Parallel, delayed
import numpy as np
import os
import math
import scipy
from pylab import *
from astropy.io import fits
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
import time
from scipy.interpolate import interp2d
from scipy.ndimage.interpolation import zoom
from scipy import optimize
from scipy import integrate
import tqdm


#os.system("taskset -p 0xff %d" % os.getpid())
# Define a series of functions to perfrom modelling of PSF
###DEPRICATED### PSF modelling moved to a different code.
# In[2]:


def gaussian(height, center_x, center_y, width_x, width_y):
     """Returns a gaussian function with the given parameters"""
     width_x = float(width_x)
     width_y = float(width_y)
     return lambda x,y: height*np.exp(
                 -(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)
 
def moments(data):
     """Returns (height, x, y, width_x, width_y)
     the gaussian parameters of a 2D distribution by calculating its
     moments """
     total = data.sum()
     X, Y = np.indices(data.shape)
     x = (X*data).sum()/total
     y = (Y*data).sum()/total
     col = data[:, int(y)]
     width_x = np.sqrt(abs((np.arange(col.size)-y)**2*col).sum()/col.sum())
     row = data[int(x), :]
     width_y = np.sqrt(abs((np.arange(row.size)-x)**2*row).sum()/row.sum())
     height = data.max()
     return height, x, y, width_x, width_y

def fitgaussian(data):
     """Returns (height, x, y, width_x, width_y)
     the gaussian parameters of a 2D distribution found by a fit"""
     params = moments(data)
     errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) -
                                  data)
     p, success = optimize.leastsq(errorfunction, params)
     return p


# Define a series of functions used to help generate non-centered galaxies

# In[3]:


def gamma(s):
    """
    Define and return the value of the gamma function.
    """
    def integrand(t, s): return t**(s-1) * np.exp(-1*t)
    gi = scipy.integrate.quad(integrand, 0, np.inf, args=s)[0]
    return gi

def NonCenSciPDF(param, xloc, yloc):
	"""
	Takes Gaussian Fit, displaces it by non integer pixel values and resamples onto original pixel grids, used for generating non-centered point Spread Funcitons 
	##DEPRICATED## PSF modelling shifted to different code.
	"""
	param2 = param	
	displacex = xloc - np.floor(xloc)
	print("disx=", displacex)
	displacey = yloc - np.floor(yloc)
	print("disy=", displacey)
	param2[1] = 29 + displacex
	param2[2] = 29 + displacey
	out = np.zeros((59,59))
	print(param)
	f = gaussian(*param2)
	for i in range (24, 34):
		for j in range (24, 34):
			out[j][i] = integrate.dblquad(f, j, j+1, lambda x: i, lambda x: i+1)[0]
	param2[1] += -displacex
	param2[2] += -displacey
	return out

def get_bn(n0):
    """
    Calculates the parameter bn from the Sersic profile.

    Args:
        n0 (int) = Sersic index.

    Returns:
        bn (float) = Value of bn.
    """
    def errfunc(bn, n):
        return abs(scipy.special.gamma(2*n0) -
                   2*scipy.special.gammainc(2*n0, bn) *
                   scipy.special.gamma(2*n0))
    bn = scipy.optimize.fmin(errfunc, 1., args=(n0,), disp=False)[0]
    return bn

def makeSersic(n0, bn0, re, ell, inc_angle, size_galaxy, disy, disx):
    """
    Calculates the flux for each pixel following a Sersic profile.

    Args:
        n0 (int) = Sersic index.
        bn0 (float) = bn parameter.
        re (float) = Effective radius in pixels.
        ell (float) = Eccentricity. Varies between 0 and 1.
        inc_angle (float) = Inclination angle in radians. Varies
                            between 0 and Pi/2.
        size_galaxy (int) = Diameter of the galaxy stamp in pixels.
    Returns:
        fl (float) = Flux for each pixel.
	disx/y = displacement of galaxy from center in pixels (between -0.5 - +0.5)
    """
    stamp = np.zeros((size_galaxy,size_galaxy))
    s2 = size_galaxy / 2
    major_axis = re
    minor_axis = re * (1.-ell)
    I_e = ((bn0)**(2*n0)) / (2*np.pi*n0*major_axis*minor_axis*gamma(2*n0))
    def f(x, y):
        x_aux = (x-s2)*np.cos(inc_angle) + (y-s2)*np.sin(inc_angle)
        y_aux = -(x-s2)*np.sin(inc_angle) + (y-s2)*np.cos(inc_angle)
        radius = np.sqrt((x_aux/major_axis)**2 + (y_aux/minor_axis)**2)
        return I_e * np.exp(-bn0*((radius)**(1./n0)))
    
    fl_tot = scipy.integrate.dblquad(f, -10000, 10000, -10000, 10000, 
                                 epsabs=1.49e-08, epsrel=1.49e-08)[0]
    print(fl_tot)
    for i in range(size_galaxy):
        def g(x):
            return i - 1./2. + disy
        def h(x):
            return i + 1./2. + disy
        for j in range(size_galaxy):
            fl = scipy.integrate.dblquad(f, j-(1./2.)+disx, j+(1./2.)+disx, g, h, 
                                 epsabs=1.49e-08, epsrel=1.49e-08)[0]
            stamp[i,j] = fl
    return stamp

def Templates(i):
	print("galaxy", i+1, "/", len(comb), comb[i])
	###GALAXY GENERATION
	t0=time.time()

	ind = comb[i][0] #Sersic Index n
	b = get_bn(ind) #Normalisation Based on n
	dx = comb[i][1]
	dy = comb[i][2]
	ell = comb[i][3]
	a = comb[i][4]
	Re= comb[i][5] * u.pc #size of galaxy half light
	Re= Re.to(u.Mpc)
	ReAng = ((Re* (1+Z)**2) / cosmo.luminosity_distance(Z)) * 57.2958 * 3600 #size of galaxy in arcsec (converting from rad) via angular distance
	#print(cosmo.luminosity_distance(Z))
	galhalflight = ReAng/res
	#print("arcsec size =", ReAng)
	#print("pixel rad =", galhalflight)
	headers = '%s' % comb[i]
	stamps = makeSersic(ind, b, galhalflight, ell, a, stampsize, -dy, -dx)
	#CountConv = CountList[i]/stamps[1]
	Galaxy = stamps 
	np.savetxt('GalaxyTemplates/GAL_Raw/%s_%s.txt' % (name, i+1), Galaxy, header=headers)
	#tot+=1
	t1=time.time()
	t01 = time.time()
	print("time taken", t1-t0)
	tottime = t01-t00
	print(i, "galaxies in", tottime, "total time") 

# Begin MAIN component of the code

# In[4]:


"""
MAIN
"""
if __name__ == '__main__':
	#tot=0
	t00 = time.time()
	### SET UP BASIC PARAMTERS AND PRIORS
	indlist = [1] #List of sersic indicies
	dispxlist = [0] # list of subpixel offsets in x and y for the center of the galaxy
	dispylist = [0] # 
	e= [0]#, 0.2, 0.4, 0.6, 0.8] #Ellipticity List
	ang = [0]#, (np.pi * 1./5.),  (np.pi * 1./5.), (np.pi * 2./5.), (np.pi * 3./5.), (np.pi * 4./5.)] #List of orientation angles
	GalSize = [3000] #np.linspace(800, 12000, 5) #half light size in parsec

	#generate a list of every combination of the above parameters
	Parameters = [indlist, dispxlist, dispylist, e, ang, GalSize]
	comb = [[]]
	for x in Parameters:
		t=[]
		for y in x:
			for i in comb:
				t.append(i+[y])
		comb = t

	print(len(comb))
	total=len(comb)
	name = 'series1_z85'

	np.savetxt('GalaxyTemplates/GAL_Raw/%s_params.txt' % name, comb)
	stampsize=101

	Z = 8.5 #Redshft

	res = 0.03 # resolution of image arcsec/pixel

	#Parallel(n_jobs=3, backend='threading', verbose=10)(delayed(Templates)(comb[i]) for i in range(total))

	with Pool(processes=3) as pool: #Parallelise work, processes is number of threads
		#pool.map(Templates, range(total))
		list(tqdm.tqdm(pool.imap(Templates, range(0,total))))

