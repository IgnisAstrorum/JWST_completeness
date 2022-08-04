#!/home/jwst/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 11:32:22 2022

@author: ignas
"""

import os
import configparser

#Read in config parameters
config = configparser.RawConfigParser()
config.read('source_insertion.config')
details_dict = dict(config.items('mask parameters'))
sextr_dict = dict(config.items('sextractor parameters'))

image_name = details_dict['img_file'] 
seg_name = details_dict['seg_file']
config_name = sextr_dict['config_file']
zeropoint = float(sextr_dict['zero_point'])

sci_img = image_name + "'[1]'"
weight_img = image_name + "'[4]'"

#Call sextractor to produce an initial segmentation file
conf_str = "sex {} -c {} -WEIGHT_IMAGE {} ".format(sci_img, config_name, weight_img)
check_string = "-CHECKIMAGE_NAME "+seg_name+" -CHECKIMAGE_TYPE SEGMENTATION"
cat_string = "-CATALOG_NAME catalog.fits -MAG_ZEROPOINT {} ".format(zeropoint)
callstring = conf_str+cat_string+check_string

os.system(callstring)

print("generating mask")
os.system("./Generate_mask.py")
print("Inserting sources")
os.system("./fake_sources.py")
print("Matching SExtractor tables")
os.system("./match_tables_data.py")
print("Outputting completeness")
os.system("./plot_completeness.py")




