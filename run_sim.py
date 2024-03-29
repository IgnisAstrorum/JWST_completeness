#!/nvme/scratch/software/anaconda3/envs/jwst/bin/python
# -*- coding: utf-8 -*-
"""
Runs the full completenes simulation
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
weight_img = image_name + "'[2]'"
if details_dict["generate_mask"]==1:
    #Call sextractor to produce an initial segmentation file 
    conf_str = "sex {} -c {} -WEIGHT_IMAGE {} ".format(sci_img, config_name, weight_img)
    check_string = "-CHECKIMAGE_NAME "+seg_name+" -CHECKIMAGE_TYPE SEGMENTATION"
    cat_string = "-CATALOG_NAME catalog.fits -MAG_ZEROPOINT {} ".format(zeropoint)
    callstring = conf_str+cat_string+check_string

    os.system(callstring)
    print("generating mask")
    os.system("./generate_mask.py")

print("Inserting sources")
os.system("./fake_sources.py")
print("Running sextractor on inserted sources")
os.system("./call_sex.py")
print("Matching SExtractor tables")
os.system("./match_tables_data.py")
print("Outputting completeness")
os.system("./plot_completeness.py")




