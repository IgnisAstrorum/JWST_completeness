# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 23:43:39 2022

@author: ignas
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

def abs_mag(app_mag, z):
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)
    lum_dist = cosmo.luminosity_distance(z).to(u.pc)
    return app_mag - 5*np.log10(lum_dist.value)+5 + 2.5*np.log10(1+z)

cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=2.725)


redshifts = np.arange(65, 175, 5)/10
print(redshifts)
comp_2d = []
for z in redshifts:
    file_temp = 'completeness_func_z{}.txt'.format(int(z*10))
    data = np.genfromtxt(file_temp, delimiter=',', skip_header=True)
    comps = []
    mags = []
    for i in data:
        comps.append(i[1])
        mags.append(abs_mag(i[0], z))

    comp_2d.append(comps)


comp_2d = np.array(comp_2d).transpose()
fig = plt.figure()
ax = fig.add_subplot(111)
contour = ax.contourf(redshifts, mags, comp_2d)
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(contour, cax=cax, orientation='vertical')
ax.set_xlabel("z")
ax.set_ylabel("Absolute magnitude")
plt.savefig('2D_completeness_z.png', dpi=300)
plt.show()




