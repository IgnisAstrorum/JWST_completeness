#!/home/jwst/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 13:36:29 2022

@author: ignas
"""
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import configparser

config = configparser.RawConfigParser()
config.read('source_insertion.config')
conf_par = dict(config.items('catalog matching'))
out_dir = conf_par['out_dir']
comp_bins = int(conf_par['comp_bins'])

comp_file = out_dir + "recovered_sources.fits"
cat_file = out_dir + "full_catalog.fits"
plot_file = out_dir + "completeness_plot.png"
func_file = out_dir + "completeness_func.txt"

completeness_sources = Table.read(comp_file)
cat_sources = Table.read(cat_file)

mag_app_comp = completeness_sources['mag_app']
mag_app_cat = cat_sources['mag_app']

det_hist, edges = np.histogram(mag_app_comp, bins=comp_bins)
cat_hist = np.histogram(mag_app_cat, bins=edges)[0]

bin_cen = [(edges[i]+edges[i+1])/2 for i in range(0, len(edges)-1)]
comp_hist = det_hist/cat_hist 

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(bin_cen, comp_hist)
ax.set_ylim(0, 1)
ax.set_xlabel("m apparent")
ax.set_ylabel("Completeness")
plt.savefig(plot_file, dpi=300)

with open(func_file, 'w') as file:
    file.write("#Completeness function mag, comp \n")
    for i in range(0, len(comp_hist)):
        string = "{},{}".format(bin_cen[i], comp_hist[i]) + "\n"
        file.write(string)







