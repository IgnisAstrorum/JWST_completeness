#!/nvme/scratch/software/anaconda3/envs/jwst/bin/python
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 10:36:17 2022
"""

import os
import numpy as np
import configparser

#Read config parameters
config = configparser.RawConfigParser()
config.read('source_insertion.config')
sex_parameters = dict(config.items('sextractor parameters'))
phys_parameters = dict(config.items('physical parameters'))
zeropoint = float(phys_parameters['zero_point'])
out_dir = sex_parameters['out_dir']
seg_fileroot = sex_parameters['seg_fileroot']
cat_fileroot = sex_parameters['cat_fileroot']
config = sex_parameters['config_file']

#Run SExtractor on the images with inserted sources
detection_files = np.genfromtxt("detection_list.txt", comments='#', dtype='str')
cat_no = 1
cat_list = "sextractor_cat.txt"
with open(cat_list, 'w') as file:
    file.write("#List of sextractor catalogs \n")
    
for file in detection_files:
    seg_file = out_dir + seg_fileroot + "{}.fits".format(cat_no)
    check_string = "-CHECKIMAGE_NAME "+seg_file+" -CHECKIMAGE_TYPE SEGMENTATION"
    in_file = file + "'[1]'"
    weights = file + "'[2]'"
    cat_name = out_dir+cat_fileroot+"{}.fits".format(cat_no)
    call1 = "sex {} -c {} -WEIGHT_IMAGE {}".format(in_file, config, weights)
    call2 = " -CATALOG_NAME {} -MAG_ZEROPOINT {} ".format(cat_name, zeropoint)
    callstring = call1 + call2+ check_string
    os.system(callstring)
    with open(cat_list, 'a') as file:
        file.write(cat_name + '\n')
    cat_no+=1



