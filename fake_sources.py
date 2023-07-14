#!/nvme/scratch/software/anaconda3/envs/jwst/bin/python
# -*- coding: utf-8 -*-
"""
Code to add sources from a stamp cataloog into a science image
"""

import numpy as np
from astropy.io import fits
import scipy.integrate as integral
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import random
from astropy.table import Table
from astropy.convolution import convolve
import configparser

def size_lum(abs_mag, ref_size=0.8, ref_mag=-21):
    return ref_size*10**(-0.2*(abs_mag-ref_mag))

def schechter (m, M = -21.19, phi = 0.0194*10**-3, a = -2.28):
    """
    Parameters
    ----------
    m : input absolute magnitude
    M : Characteristic magnitude The default is -21.19.
    phi : Scaling The default is 0.0194*10**-3.
    a : slope The default is -2.28.
    
    Calculates the value of the schechter LF at magnitude m
    """
    return phi*np.log(10)/2.5*10**(-0.4*(m - M)*(a+1))*np.exp(-10**(-0.4*(m - M)))

def dist_modulus(z, cosmo):
    """
    Parameters
    ----------
    z : Redshift
    cosmo : Parameters of the cosmology used
    
    Calculates the distance modulus at redshift z for a particular cosmology
    """
    lum = cosmo.luminosity_distance(z)/u.Mpc*10**6
    return 5*np.log10(lum)-5 - 2.5*np.log10(1+z)

def random_magnitude(mag_max, mag_min, n_sample, use_LF, lf_params):
    """
    Parameters
    ----------
    mag_max : maximum magnitude that can be picked
    mag_min : minimum magnitude that can be picked
    n_sample : # of magnitude values that can be chosen from
    Returns
    -------
    A random absolute magnitude from a chosen range.
    """
    mag_range = np.linspace(mag_max, mag_min, n_sample)
    bin_width = (mag_min - mag_max)/n_sample
    if use_LF:
        mag_char = float(lf_params['mag_char'])
        scale = float(lf_params['scale'])
        slope = float(lf_params['slope'])
        args = (mag_char, scale, slope)
        norm = integral.quad(schechter, mag_max, mag_min, args = args)[0] #Normalize LF
        #Generate list of magnitude probabilities
        mag_probs = schechter(mag_range)*bin_width/norm
        mag_probs /= sum(mag_probs) #Adjust probabilities for numpy
        return np.random.choice(mag_range, p=mag_probs)
    else: return np.random.choice(mag_range)

def get_fits_data(filename):
    """
    Parameters
    ----------
    filename : Name of file to open (string)
    Returns
    -------
    data : array containing the read-in data
    """
    hdul = fits.open(filename)
    data = hdul[0].data
    hdul.close()
    return data

def open_stamp_params(param_file):
    """
    Parameters
    ----------
    param_file : A string containing the location of the stamp parameter file
    Returns
    -------
    param_array : Array containing stamp parameters for each file #
    """
    file_params = []
    with open(param_file, "r") as params:
        number = 1
        for line in params.readlines():
            file_vals = [float(i) for i in line.strip().split(" ")]
            file_vals.append(number)
            number+=1
            file_params.append(file_vals)
    
    param_array = np.array(file_params)
    return param_array

def generate_catalog(mag_max, mag_min, phys_params, lf_params, 
                   file_params, cosmology, cat_list):
    """
    Parameters
    ----------
    segmentation_file : str
        Name of the segmentation file
    image_file : str
        Name of the image file

    Returns
    -------
    Generates a list of sources with different mag_app and coordinates for insertion

    """
    #read in config parameters
    dir_path = file_params['stamp_dir']
    mask_name = file_params['mask_file']
    n_inserts = int(phys_params['n_inserts'])
    use_LF = int(phys_params['use_lf'])
    H0 = float(cosmology['h0'])
    omega_m = float(cosmology['omega_m'])
    sizes = np.array(phys_params["galaxy_sizes"].split(',')).astype(float)
    zero_point = float(phys_params['zero_point'])
    outroot = file_params['outroot']
    stamp_params = file_params['stamp_params']
    z = float(phys_params['redshift'])
    r_ref = float(phys_params['r_ref'])
    m_ref = float(phys_params['m_ref'])
    
    #Open catalog parameters for file choice
    param_file = dir_path + stamp_params 
    param_array = open_stamp_params(param_file)
    
    #Set lambdacdm cosmology        
    cosmo = FlatLambdaCDM(H0, omega_m)
    
    #Get mask data
    mask = get_fits_data(mask_name)

    free_coords = np.where(mask==0)
    
    y_vals = free_coords[0] 
    x_vals = free_coords[1]
    
    #Generate source catalog
    source_catalog = []
    #Coords
    coordinates = [[],[]]
    for i in range (0, n_inserts):
        #Generate random unoccupied position
        position_index = random.randint(0, len(x_vals)-1)
        x = x_vals[position_index]
        y = y_vals[position_index]
        """
        while x > 5000:
            position_index = random.randint(0, len(x_vals)-1)
            x = x_vals[position_index]
            y = y_vals[position_index]"""
        
        #Sample magnitudes, convert to apparent
        abs_mag = random_magnitude(mag_max, mag_min, 20000, use_LF, lf_params)
        modulus = float(dist_modulus(z, cosmo))
        vis_mag = abs_mag+modulus
        tot_counts = 10**(-(vis_mag - zero_point)/2.5)
        
        #Generate stamp size
        #Use a size-luminosity relation
        size_mag = size_lum(abs_mag, ref_size = r_ref, ref_mag = m_ref)
        size_diff = np.abs(sizes - size_mag)
        size = sizes[np.where(size_diff == min(size_diff))[0]]*1000
        
        #Pick file according to size
        size_rows = np.where(param_array[:,-2] == size)[0]
        file_nos = [i[-1] for i in param_array[size_rows]]
        file_no = int(np.random.choice(file_nos))
        
        #Add source to the catalog
        source_catalog.append([int(x), int(y), z, vis_mag, tot_counts, file_no, size, abs_mag])
        coordinates[0].append(x)
        coordinates[1].append(y)

    #Save catalog
    table = Table([[i[0] for i in source_catalog], [i[1] for i in source_catalog], 
                   [i[2] for i in source_catalog], [i[3] for i in source_catalog], 
                   [i[4] for i in source_catalog], [i[5] for i in source_catalog],
                  [i[6] for i in source_catalog], [i[7] for i in source_catalog]],
                  names=('x', 'y', 'redshift', 'mag_app', 'counts', 'file_no', 
                         'size', 'abs_mag'))
    centre_bin = -(mag_max + mag_min)/2
    catalog_filename = outroot+'catalog_{0:.1f}.fits'.format(centre_bin)
    table.write(catalog_filename, format='fits', overwrite=True)
    with open(cat_list, 'a') as file:
        file.write(catalog_filename + '\n')
    return source_catalog


def insert_sources(mag_max, mag_min, phys_params, lf_params, 
                   file_params, cosmology, detection_list, cat_list):
    """
    Parameters
    ----------
    mag_max : maximum magnitude
    mag_min : minimum magnitude
    phys_params : dictionary of physical galaxy parameters
    lf_params : parameters for Schechter LF
    file_params : names for i/o files
    cosmology : parameters of the lambdaCDM
    detection_list : filename for file containing the sextractor image list
    
    Generates a catalog of sources between specified mag_min and mag_max
    and inserts the sources into an image file

    """
    source_catalog = generate_catalog(mag_max, mag_min, phys_params, lf_params, 
                       file_params, cosmology, cat_list)
    
    image_file = file_params['image_file']
    psf_file = file_params['psf_file']
    stamp_size = int(phys_params['stamp_size'])
    stamp_dir = file_params['stamp_dir']
    outroot = file_params['outroot']
    fileroot = file_params['stamp_fileroot']
    
    #Open simulated data file
    image = fits.open(image_file, memmap=False)
    image_data = image['SCI'].data
    #Set PSF as convolution kernel
    psf_file = fits.open(psf_file)
    kernel = psf_file[0].data
    #Trim kernel down to stamp size
    center = kernel.shape[0]//2
    kernel = kernel[center-stamp_size//2:center+stamp_size//2+1, center-stamp_size//2:center+stamp_size//2+1]
    psf_file.close()



    #Add sources
    radius = int((stamp_size-1)/2) #stamp radius
    for source in source_catalog:
        file_no = source[5]
        stamp_file = stamp_dir+fileroot+"{}.txt".format(file_no)
        stamp_data = np.genfromtxt(stamp_file, skip_header=1)
        stamp_norm = stamp_data/np.sum(stamp_data)*source[4]
        stamp_corr = convolve(stamp_norm, kernel)
        image_data[source[1]-radius:source[1]+radius+1, 
               source[0]-radius:source[0]+radius+1] += stamp_corr
    
    centre_bin2 = -(mag_max + mag_min)/2 #Write changes
    with open(detection_list, 'a') as file:
        file.write(outroot+"added_{0:.1f}.fits".format(centre_bin2) + '\n')
    image.writeto(outroot+"added_{0:.1f}.fits".format(centre_bin2), overwrite=True)
    image.close()


def main():
    config = configparser.RawConfigParser()
    config.read('source_insertion.config')
    
    #Read in galaxy parameters from config
    phys_params = dict(config.items('physical parameters'))
    mag_min = float(phys_params['mag_min'])
    mag_max = float(phys_params['mag_max'])
    bin_width = float(phys_params['bin_width'])
    #Read in LF parameters from config
    lf_params = dict(config.items('LF parameters'))
    #Read in output file parameters
    file_params = dict(config.items('files'))
    #Read in cosmology params
    cosmology = dict(config.items('cosmology'))
    
    bin_edges = [mag_max + i*bin_width for i in range(0, int((mag_min - mag_max)/bin_width))]
    
    file_list = "detection_list.txt"
    with open(file_list, 'w') as file:
        file.write('#List of sextractor files \n')
        
    cat_list = "generated_cat.txt"
    with open(cat_list, 'w') as file:
        file.write('#List of generated cat files \n')
    for i in range(0, len(bin_edges)-1):  
        insert_sources(bin_edges[i], bin_edges[i+1], phys_params,
                       lf_params, file_params, cosmology, file_list, cat_list)

if __name__ == '__main__':
    main()

