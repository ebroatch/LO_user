"""
Test plot to get subplot spacing correct

"""
import numpy as np
import xarray as xr
import pickle
from datetime import datetime, timedelta
import pandas as pd
from cmocean import cm

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
#import pinfo
from importlib import reload
reload(pfun)
#reload(pinfo)

Ldir = Lfun.Lstart()
if '_mac' in Ldir['lo_env']: # mac version
    pass
else: # remote linux version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt

# START
#ds = xr.open_dataset(in_dict['fn'])
aa1 = [-2, 0, 44.5, 46.5]
aa2 = [-4, 4, 43, 47]
aa3 = [-0.2, 1.2, 44.9, 45.1]
# find aspect ratio of the map
# AR is the aspect ratio of the map: Vertical/Horizontal
# AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
# fs = 14
# hgt = 10
#pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
x = np.linspace(-4,4)
y = np.linspace(43,47)
XX, YY = np.meshgrid(x, y)
vort = 2*XX+3*(YY-45)

# set color limits
vv = 2*np.nanstd(vort)
vmin = -vv
vmax = vv

fig = plt.figure(figsize=(14,8))
#gs = fig.add_gridspec(nrows=4, ncols=5, hspace=0.5, wspace=0.4)
gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[17,10,1], height_ratios=[3,1])
cmap = cm.curl

ax1 = fig.add_subplot(gs[0,1]) 
cs1 = plt.pcolormesh(XX, YY, vort, cmap=cmap, vmin = vmin, vmax = vmax)
ax1.set_title('Plume focus', fontsize=12)
#fig.colorbar(cs1)
ax1.axis(aa1)
pfun.dar(ax1)
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')

ax2 = fig.add_subplot(gs[0,0])
cs2 = plt.pcolormesh(XX, YY, vort, cmap=cmap, vmin = vmin, vmax = vmax)
ax2.set_title('Full model', fontsize=12)
#fig.colorbar(cs2)
ax2.axis(aa2)
pfun.dar(ax2)
ax2.set_xlabel('Longitude')
ax2.set_ylabel('Latitude')
#pfun.add_info(ax2, in_dict['fn'])
# pfun.add_coast(ax2)
# pfun.add_bathy_contours(ax, ds, txt=True)

ax3 = fig.add_subplot(gs[1,0:2])
cs3 = plt.pcolormesh(XX, YY, vort, cmap=cmap, vmin = vmin, vmax = vmax)
ax3.set_title('Estuary focus', fontsize=12)
ax3.axis(aa3)
pfun.dar(ax3)
ax3.set_xlabel('Longitude')
ax3.set_ylabel('Latitude')

ax4 = fig.add_subplot(gs[:,2])
plt.colorbar(cs3, cax=ax4)
plt.suptitle('Surface Vorticity $[s^{-1}]$', fontsize=16)
#plt.tight_layout()

plt.show()
