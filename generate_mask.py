#!/nvme/scratch/software/anaconda3/envs/jwst/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:19:49 2022

@author: Ignas
A simple program to generate a mask from a segmentation map and the original science image
This program requires the configuration file mask.config to run.
Parameters to be specified:
    SExtractor segmentation file name
    Science image file name
    Edge buffer size (stamp size if used for completeness simulations)
    Object buffer size
    Output file name
"""

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
import configparser

def load_config():
    config = configparser.RawConfigParser()
    config.read('source_insertion.config')
    details_dict = dict(config.items('mask parameters'))
    return details_dict
    
def main():
    #Load configuration parameters
    generate_mask, segmentation_file, image, stamp_siz, object_buffer, output_filename = load_config()
    par_dict = load_config()
    segmentation_file = par_dict['seg_file']
    image = par_dict['img_file']
    stamp_siz = int(par_dict['stamp_size'])
    object_buffer = int(par_dict['buffer'])
    output_filename = par_dict['out_file']

    #Open the image and the segmentation file
    imag = fits.open(image, memmap=False)
    im_data = imag[1].data
    imag.close()

    seg = fits.open(segmentation_file, memmap=False)
    seg_data = seg[0].data
    seg.close()

    mask = np.zeros((len(im_data), len(im_data[0])))


    #mask blank areas and a buffer around them
    x0, y0 = np.where(im_data == 0)
    mask[x0, y0] = -2
    maskedx, maskedy = np.where(mask == -2)
    allx, ally = np.where(mask==mask)
    #Use astropy catalog matching to determine pixel distances
    mask_xy = SkyCoord(maskedx, maskedy, np.zeros(len(maskedy)),
                   representation_type='cartesian', unit='pix')
    all_xy = SkyCoord(allx, ally, np.zeros(len(ally)), 
                      representation_type='cartesian', unit='pix')
    match_id, d2d, d3d = all_xy.match_to_catalog_3d(mask_xy)
    criterion = d3d<(100)*u.pix
    #Mask the buffer zone
    extra_mask = all_xy[criterion]
    extra_x = extra_mask.x.value.astype(np.int64)
    extra_y = extra_mask.y.value.astype(np.int64)
    mask[extra_x, extra_y] = 1
    
    #Mask objects
    xseg, yseg = np.where(seg_data != 0)
    mask[xseg, yseg] = -3
    #Mask a buffer zone around objects if desired
    if object_buffer > 0:
        objy, objx = np.where(mask==-3)
        obj_xy = mask_xy = SkyCoord(objy, objx, np.zeros(len(objy)),
                       representation_type='cartesian', unit='pix')
        obj_id, obj_2d, obj_3d = all_xy.match_to_catalog_3d(obj_xy)
        buffer_crit = obj_3d<object_buffer*u.pix
        buff_mask = all_xy[buffer_crit]
        buff_x = buff_mask.x.value.astype(np.int64)
        buff_y = buff_mask.y.value.astype(np.int64)
        mask[buff_x, buff_y] = 1
        
    #mask image edges
    mask+=1
    mask[stamp_siz:len(mask)-stamp_siz-1, stamp_siz:len(mask[0])-stamp_siz-1] -= 1
    mask[np.where(mask!=0)[0], np.where(mask!=0)[1]] = 1
    plt.imshow(mask)

    mask_hdul = fits.PrimaryHDU(mask)
    mask_hdul.writeto(output_filename, overwrite=True)
    
if __name__ == '__main__':
    main()


