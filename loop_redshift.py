#!/nvme/scratch/software/anaconda3/envs/jwst/bin/python
# -*- coding: utf-8 -*-
"""
A basic code to loop the simulation across a redshift range
"""

import os
import numpy as np
import configparser

z_list = np.arange(6.5, 17.5, 0.5)
print(z_list)

for z in z_list:
    config = configparser.ConfigParser()
    config.read('source_insertion.config')
    phys_par = dict(config.items("physical parameters"))
    phys_par["redshift"] = z
    config["physical parameters"] = phys_par
    with open("source_insertion.config", 'w') as conf:
        config.write(conf)
    
    os.system("./run_sim.py")
