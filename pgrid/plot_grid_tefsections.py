"""
Plot grid to have a look at it. Accepts an optional command line argument
to look at a grid other than the one set in gfun.py.
"""
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pickle
from lo_tools import Lfun
from lo_user_tools import llxyfun #my functions for converting lat/lon to m

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gridname', default='',type=str) # e.g. cas6
parser.add_argument('-dmax', default=5, type=int) # max depth for colormap [m]
parser.add_argument('-small', default=False, type=Lfun.boolean_string) # True for laptop size
args = parser.parse_args()
zmin = -args.dmax

import gfun
if len(args.gridname) > 0:
    Gr = gfun.gstart(gridname=args.gridname)
else:
    Gr = gfun.gstart()
from lo_tools import plotting_functions as pfun
import gfun_plotting as gfp


testing = True
if testing:
    from importlib import reload
    reload(gfun)
    reload(gfp)

# select grid file
in_fn = gfun.select_file(Gr)

# load the default choices
try:
    dch = pickle.load(open(Gr['gdir'] / 'choices.p', 'rb'))
except FileNotFoundError:
    # you could fill this in by hand if you wanted
    dch = {'analytical': False} # hack to make cas6 work

# get river info if it exists
do_riv = False
ri_fn = Gr['gdir'] / 'roms_river_info.csv'
if ri_fn.is_file():
    rri_df = pd.read_csv(ri_fn, index_col='rname')
    do_riv = True

# load the data
ds = xr.open_dataset(in_fn)
z = -ds.h.values
mask_rho = ds.mask_rho.values

lon = ds.lon_rho.values
lat = ds.lat_rho.values

plon, plat = pfun.get_plon_plat(lon,lat)
pad = 0.05*(plat[-1,0]-plat[0,0])
ax_lims = (plon[0,0]-pad, plon[0,-1]+pad, plat[0,0]-pad, plat[-1,0]+pad)

# make a version of z with nans where masked
zm = z.copy()
zm[mask_rho == 0] = np.nan

# PLOTTING
plt.close('all')
if args.small:
    figsize = (8,8)
else:
    figsize = (12,12)
pfun.start_plot(figsize=figsize)

# bathymetry
fig = plt.figure()
ax = fig.add_subplot(111)
#cs = ax.pcolormesh(llxyfun.lon2x(lon,0,45), llxyfun.lat2y(lat,45), zm, vmin=-200, vmax=0, cmap='Spectral_r') #change to plot on distance grid instead of lat-lon and adjust color limits
cs = ax.pcolormesh(plon, plat, zm, vmin=-300, vmax=0, cmap='Spectral_r')
# cs = ax.pcolormesh(plon, plat, zm, vmin=-120, vmax=-100, cmap='Spectral_r')
#fig.colorbar(cs, ax=ax)
if dch['analytical'] == True:
    pass
else:
    pfun.add_coast(ax)
pfun.dar(ax)
#ax.axis(ax_lims)
ax.set_xlim(-0.2,2.5)
ax.set_ylim(44.9,45.1)
#ax.set_title(in_fn.name)
#ax.set_title('Idealized model bathymetry with 20km sill')
if args.gridname=='sill20kmdeep':
    ax.set_title('20km sill')
    sect_lon=llxyfun.x2lon([40,41.25,42.5,43.75,45],0,45)
elif args.gridname=='sill5km':
    ax.set_title('5km sill')
    sect_lon=llxyfun.x2lon([40,45,50,55,60],0,45)
elif args.gridname=='sill80km':
    ax.set_title('80km sill')
    sect_lon=llxyfun.x2lon([40,60,80,100,120],0,45)

for i in range(5):
    ax.plot([sect_lon[i],sect_lon[i]],[44.95,45.05],'k',lw=2)
#ax.text(.05, .95, Gr['gridname'], transform=ax.transAxes)
#ax.text(.95, .05, str(mask_rho.shape), ha='right', transform=ax.transAxes, bbox=pfun.bbox)
if do_riv:
    gfp.add_river_tracks(Gr, ds, ax)

if False:    
    # mask
    fig = plt.figure()
    ax = fig.add_subplot(111)
    tt = ['rho', 'u', 'v']
    sym = dict(zip(['rho', 'u', 'v'],['o','>','^']))
    c = dict(zip(['rho', 'u', 'v'],['b','orange','r']))
    for t in tt:
        x = ds['lon_'+t].values
        y = ds['lat_'+t].values
        m = ds['mask_'+t].values
        ax.plot(x, y, sym[t], c=c[t], alpha=.2, ms=3)
        ax.plot(x[m==1], y[m==1], sym[t], c=c[t], ms=3)
    
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(ax_lims)
    ax.set_title(in_fn.name)
    ax.text(.05, .95, Gr['gridname'], transform=ax.transAxes)
        
ds.close()

plt.show()
