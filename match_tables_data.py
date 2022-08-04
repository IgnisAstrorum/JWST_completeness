#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Matches SExtractor catalogs to inserted catalogs
"""

import numpy as np
from astropy.table import Table, join, Column, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
import configparser

config = configparser.RawConfigParser()
config.read('source_insertion.config')
conf_par = dict(config.items("catalog matching"))
error = float(conf_par['error'])
out_dir = conf_par['out_dir']

#Output table parameters
comp_tables = []
comp_tablename = out_dir + "recovered_sources.fits"

cat_tables = []
cat_tablename = out_dir + "full_catalog.fits"
#Table for undetected inserted sources
ins_und = []
und_tablename = out_dir + "inserted_undetected.fits"
sextractor_files = np.genfromtxt("sextractor_cat.txt", comments='#', dtype='str')
inserted_cat_files = np.genfromtxt("generated_cat.txt", comments='#', dtype='str')
for index in range(0, len(sextractor_files)):
    #Get cartesian and eq coordinates of recovered and generated catalogs
    extracted_file = sextractor_files[index]

    data_ext = Table.read(extracted_file)# extracted_sources[1].data
    extracted_xy = SkyCoord(data_ext['X_IMAGE'], data_ext['Y_IMAGE'],
                            np.zeros(len(data_ext['X_IMAGE'])),
                            representation_type='cartesian', unit='pix')
    extracted_sky = SkyCoord(ra=data_ext['ALPHA_J2000'],
                             dec=data_ext['DELTA_J2000'], unit='deg')

    cat_file = inserted_cat_files[index]
    cat_data = Table.read(cat_file)#cat_sources[1].data
    cat_tables.append(cat_data)
    catalog_xy = SkyCoord(cat_data['x'], cat_data['y'],
                            np.zeros(len(cat_data['y'])),
                            representation_type='cartesian', unit='pix')

    #Match recovered sources to fake catalog
    match_id, d2d, d3d = extracted_xy.match_to_catalog_3d(catalog_xy)
    constraint = d3d <= error*u.pix
    
    dist_col = Column(d3d[constraint], name='NN distance')
    #Join by unique ID to avoid extra rows
    ID_col = Column([i for i in range(0, len(d3d[constraint]))], name='ID')
    matches_cat = cat_data[match_id[constraint]]
    matches_cat.add_column(dist_col)
    matches_cat.add_column(ID_col)
    ext_matches = data_ext[constraint]
    ext_matches.add_column(ID_col)
    
    #print(len(matches_cat))
    if len(matches_cat) > 0:
        recovered_inserted = join(ext_matches, matches_cat, keys='ID')
        comp_tables.append(recovered_inserted)
    
    und_id, und_dist2, und_dist3 = catalog_xy.match_to_catalog_3d(extracted_xy)
    und_cond = und_dist3 > error*u.pix
    und_cat = cat_data[und_cond]
    ins_und.append(und_cat)
    


#stack recovered tables
final_und_table = ins_und[-1]
und_co = SkyCoord(final_und_table['x'], final_und_table['y'],
                        np.zeros(len(final_und_table['y'])),
                        representation_type='cartesian', unit='pix')
for i in range(1, len(ins_und)):
    match_co_u = SkyCoord(ins_und[i]['x'], ins_und[i]['y'],
                            np.zeros(len(ins_und[i]['x'])),
                            representation_type='cartesian', unit='pix')
    print(len(und_co))
    matched, d2d_u, d3d_u = match_co_u.match_to_catalog_3d(und_co)
    new_cond_u = d3d_u > 0*u.pix
    final_und_table = vstack([final_und_table, ins_und[i][new_cond_u]])

final_und_table.write(und_tablename, format='fits', overwrite=True)


final_comp_table = comp_tables[0]
comp_co = SkyCoord(final_comp_table['X_IMAGE'], final_comp_table['Y_IMAGE'],
                        np.zeros(len(final_comp_table['X_IMAGE'])),
                        representation_type='cartesian', unit='pix')
for i in range(1, len(comp_tables)):
    match_co_c = SkyCoord(comp_tables[i]['X_IMAGE'], comp_tables[i]['Y_IMAGE'],
                            np.zeros(len(comp_tables[i]['X_IMAGE'])),
                            representation_type='cartesian', unit='pix')
    matched, d2d_c, d3d_c = match_co_c.match_to_catalog_3d(comp_co)
    new_cond_c = d3d_c > 0*u.pix
    final_comp_table = vstack([final_comp_table, comp_tables[i][new_cond_c]])
final_comp_table.write(comp_tablename, format='fits', overwrite=True)


full_catalog = cat_tables[0]
for i in range(1, len(cat_tables)):
    full_catalog = vstack([full_catalog, cat_tables[i]])
full_catalog.write(cat_tablename, format='fits', overwrite=True)



