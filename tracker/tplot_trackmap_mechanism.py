"""
Plot results of a particle tracking experiment.
"""

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
from lo_user_tools import llxyfun
Ldir = Lfun.Lstart()

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

# Choose an experiment and release to plot.
in_dir0 = Ldir['LOo'] / 'tracks'
exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    itext='** Choose experiment from list **', last=False)
rel = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)

# get Datasets
fn = in_dir0 / exp_name / rel
fng = in_dir0 / exp_name / 'grid.nc'
dsr = xr.open_dataset(fn, decode_times=False)
dsg = xr.open_dataset(fng)

NT, NP = dsr.lon.shape

# get a list of datetimes
ot_vec = dsr.ot.values
dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

# gather some fields, for convenience
lonp, latp = pfun.get_plon_plat(dsg.lon_rho.values, dsg.lat_rho.values)
hh = dsg.h.values
maskr = dsg.mask_rho.values
#

#sillmid = llxyfun.x2lon(44e3,0,45)
# sillsea = llxyfun.x2lon(40e3,0,45)
# sillland = llxyfun.x2lon(60e3,0,45)

#grid='20kmdeep'
sillsea = llxyfun.x2lon(40e3,0,45)
#if grid=='5km':
if exp_name.split('_')[0]=='sill5kmest':
    sillland = llxyfun.x2lon(45e3,0,45)
    #xlonlim=1.1
    aa=[0,1.1,44.95,45.05] #estuary focus limits
#elif grid=='20kmdeep':
elif exp_name.split('_')[0]=='sill20kmdeepest':
    sillland = llxyfun.x2lon(60e3,0,45)
    #xlonlim=1.3
    aa=[0,1.3,44.95,45.05] #estuary focus limits
# elif grid=='80km':
elif exp_name.split('_')[0]=='sill80kmest':
    sillland = llxyfun.x2lon(120e3,0,45)
    #xlonlim=2.1
    aa=[0,2.1,44.95,45.05] #estuary focus limits
lon1 = dsr.lon.where((dsr.lon.sel(Time=0)<sillland) & (dsr.lon.sel(Time=0)>sillsea),drop=True).values
lat1 = dsr.lat.where((dsr.lon.sel(Time=0)<sillland) & (dsr.lon.sel(Time=0)>sillsea),drop=True).values
lon2 = dsr.lon.where((dsr.lon.sel(Time=12)<sillland) & (dsr.lon.sel(Time=12)>sillsea),drop=True).values
lat2 = dsr.lat.where((dsr.lon.sel(Time=12)<sillland) & (dsr.lon.sel(Time=12)>sillsea),drop=True).values

# lon = dsr.lon.values[:,::step]
# lat = dsr.lat.values[:,::step]
# lon = dsr.lon.values #use this for no subsampling
# lat = dsr.lat.values

# make a mask that is False from the time a particle first leaves the domain
# and onwards (outer basin particles)
AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
        dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
ib_mask = np.ones(lon1.shape, dtype=bool)
ib_mask[lon1 < AA[0]] = False
ib_mask[lon1 > AA[1]] = False
ib_mask[lat1 < AA[2]] = False
ib_mask[lat1 > AA[3]] = False
NTS, NPS = lon1.shape
for pp in range(NPS):
    tt = np.argwhere(ib_mask[:,pp]==False)
    if len(tt) > 0:
        ib_mask[tt[0][0]:, pp] = False

# and apply the mask to lon and lat
lon1[~ib_mask] = np.nan
lat1[~ib_mask] = np.nan

# make a mask that is False from the time a particle first leaves the domain 
# and onwards (inner basin particles)
AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
        dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
ib_mask = np.ones(lon2.shape, dtype=bool)
ib_mask[lon2 < AA[0]] = False
ib_mask[lon2 > AA[1]] = False
ib_mask[lat2 < AA[2]] = False
ib_mask[lat2 > AA[3]] = False
NTS, NPS = lon2.shape
for pp in range(NPS):
    tt = np.argwhere(ib_mask[:,pp]==False)
    if len(tt) > 0:
        ib_mask[tt[0][0]:, pp] = False

# and apply the mask to lon and lat
lon2[~ib_mask] = np.nan
lat2[~ib_mask] = np.nan

# PLOTTING - SPAGHETTI PLOT
plt.close('all')
pfun.start_plot(figsize=(10,10))
# fig = plt.figure()
fig, [ax,ax2] = plt.subplots(2,1,figsize=(20,10))

# MAP
# set domain limits
# if False:
#     # plot full domain
#     aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
# else:
#     # automatically plot region of particles, with padding #use outer particles to pick limits
#     pad = .02
#     aa = [np.nanmin(lon1) - pad, np.nanmax(lon1) + pad,
#     np.nanmin(lat1) - pad, np.nanmax(lat1) + pad]
# aa=[0,1.3,44.95,45.05] #estuary focus limits
    
#ax = fig.add_subplot(121)
#ax = fig.add_subplot(111)
zm = -np.ma.masked_where(maskr==0, hh)
ax.pcolormesh(lonp, latp, zm, vmin=-300, vmax=20,
    cmap='Greys_r')
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title('12h tracks starting on sill')
ax2.pcolormesh(lonp, latp, zm, vmin=-300, vmax=20,
    cmap='Greys_r') #change vmin to -300 to cut out blackest part of cmap
ax2.axis(aa)
pfun.dar(ax2)
ax2.set_xlabel('Longitude')
ax2.set_ylabel('Latitude')
ax2.set_title('12h tracks ending on sill')
# add the tracks (packed [time, particle])
# regular spaghetti plots

# # subsample output for plotting #SKIP SUBSAMPLING
# npmax = 600 # max number of points to plot
# step = max(1,int(np.floor(NP/npmax))) #use same for both assuming approx equal number of particles in each basin
step=1 #no subsampling
ax.plot(lon1[:13,::step], lat1[:13,::step], '-', color='tab:cyan', linewidth=.2) #DON'T PLOT LINES FOR NOW #comment out to skip lines
ax2.plot(lon2[:13,::step], lat2[:13,::step], '-', color='tab:pink', linewidth=.2)
# ax.plot(lon[0,:], lat[0,:], 'og', alpha=.3)
# ax.plot(lon[-1,:], lat[-1,:], 'or', alpha=.3)
#ax.plot(lon[0,:], lat[0,:], '.g', alpha=.3, markeredgecolor='none')
#ax.plot(lon[-1,:], lat[-1,:], '.r', alpha=.3, markeredgecolor='none')
ax.plot(lon1[12,:], lat1[12,:], '.', color='tab:blue', alpha=1, markeredgecolor='none')
ax2.plot(lon2[0,:], lat2[0,:], '.', color='tab:red', alpha=1, markeredgecolor='none')

# # time series
# td = (ot_vec - ot_vec[0])/86400
# tv_list = ['z', 'u', 'v']
# #tv_list = ['u', 'v', 'lon', 'lat']
# ntv = len(tv_list)
# for ii in range(ntv):
#     tv = tv_list[ii]
#     NC = 2
#     ax = fig.add_subplot(ntv,NC, (ii+1)*NC)
#     # v = dsr[tv].values[:,::step]
#     v = dsr[tv].values
#     v[~ib_mask] = np.nan
#     ax.plot(td, v, lw=.5, alpha=.2)
#     ax.text(.05, .05, tv, fontweight='bold', transform=ax.transAxes)
#     if ii == ntv-1:
#         ax.set_xlabel('Time (days)')

#plt.show()
fn_fig = Ldir['LOo'] / 'plots' / 'tplot_trackmap_mechanism.png'
plt.savefig(fn_fig)
#plt.show()
pfun.end_plot()


#NEW PLOT DIFFERENT SORTING
fig, [ax0,ax1,ax2,ax3] = plt.subplots(4,1,figsize=(20,20))

#particles that start and end on the sill
lon0 = dsr.lon.where((dsr.lon.sel(Time=0)<sillland) & (dsr.lon.sel(Time=0)>sillsea) & (dsr.lon.sel(Time=12)<sillland) & (dsr.lon.sel(Time=12)>sillsea),drop=True).values
lat0 = dsr.lat.where((dsr.lon.sel(Time=0)<sillland) & (dsr.lon.sel(Time=0)>sillsea) & (dsr.lon.sel(Time=12)<sillland) & (dsr.lon.sel(Time=12)>sillsea),drop=True).values
#particles that start on sill and end off sill
lon1 = dsr.lon.where((dsr.lon.sel(Time=0)<sillland) & (dsr.lon.sel(Time=0)>sillsea) & ( (dsr.lon.sel(Time=12)>sillland) | (dsr.lon.sel(Time=12)<sillsea) ),drop=True).values
lat1 = dsr.lat.where((dsr.lon.sel(Time=0)<sillland) & (dsr.lon.sel(Time=0)>sillsea) & ( (dsr.lon.sel(Time=12)>sillland) | (dsr.lon.sel(Time=12)<sillsea) ),drop=True).values
#particles that start off sill and end on sill
lon2 = dsr.lon.where(( (dsr.lon.sel(Time=0)>sillland) | (dsr.lon.sel(Time=0)<sillsea) ) & (dsr.lon.sel(Time=12)<sillland) & (dsr.lon.sel(Time=12)>sillsea),drop=True).values
lat2 = dsr.lat.where(( (dsr.lon.sel(Time=0)>sillland) | (dsr.lon.sel(Time=0)<sillsea) ) & (dsr.lon.sel(Time=12)<sillland) & (dsr.lon.sel(Time=12)>sillsea),drop=True).values
#particles that switch basins
lon3 = dsr.lon.where(( (dsr.lon.sel(Time=0)>sillland) & (dsr.lon.sel(Time=12)<sillsea) ) | (  (dsr.lon.sel(Time=0)<sillsea) & (dsr.lon.sel(Time=12)>sillland) ),drop=True).values
lat3 = dsr.lat.where(( (dsr.lon.sel(Time=0)>sillland) & (dsr.lon.sel(Time=12)<sillsea) ) | (  (dsr.lon.sel(Time=0)<sillsea) & (dsr.lon.sel(Time=12)>sillland) ),drop=True).values

zm = -np.ma.masked_where(maskr==0, hh)
ax0.pcolormesh(lonp, latp, zm, vmin=-300, vmax=20, cmap='Greys_r')
ax1.pcolormesh(lonp, latp, zm, vmin=-300, vmax=20, cmap='Greys_r')
ax2.pcolormesh(lonp, latp, zm, vmin=-300, vmax=20, cmap='Greys_r')
ax3.pcolormesh(lonp, latp, zm, vmin=-300, vmax=20, cmap='Greys_r')

ax0.set_xlabel('Longitude')
ax0.set_ylabel('Latitude')
ax1.set_xlabel('Longitude')
ax1.set_ylabel('Latitude')
ax2.set_xlabel('Longitude')
ax2.set_ylabel('Latitude')
ax3.set_xlabel('Longitude')
ax3.set_ylabel('Latitude')
ax0.set_title('12h tracks starting and ending on sill')
ax1.set_title('12h tracks starting on sill and ending off sill')
ax2.set_title('12h tracks starting off sill and ending on sill')
ax3.set_title('12h tracks switching basins')

ax0.plot(lon0[:13,::step], lat0[:13,::step], '-', color='tab:blue', linewidth=.1, label='Track') #plot hours 0 to 12
ax0.plot(lon0[0,:], lat0[0,:], '.', color='tab:green', label='Start')
ax0.plot(lon0[12,:], lat0[12,:], '.', color='tab:red', label='End')
ax1.plot(lon1[:13,::step], lat1[:13,::step], '-', color='tab:blue', linewidth=.1, label='Track') #plot hours 0 to 12
ax1.plot(lon1[0,:], lat1[0,:], '.', color='tab:green', label='Start')
ax1.plot(lon1[12,:], lat1[12,:], '.', color='tab:red', label='End')
ax2.plot(lon2[:13,::step], lat2[:13,::step], '-', color='tab:blue', linewidth=.1, label='Track') #plot hours 0 to 12
ax2.plot(lon2[0,:], lat2[0,:], '.', color='tab:green', label='Start')
ax2.plot(lon2[12,:], lat2[12,:], '.', color='tab:red', label='End')
ax3.plot(lon3[:13,::step], lat3[:13,::step], '-', color='tab:blue', linewidth=.1, label='Track') #plot hours 0 to 12
ax3.plot(lon3[0,:], lat3[0,:], '.', color='tab:green', label='Start')
ax3.plot(lon3[12,:], lat3[12,:], '.', color='tab:red', label='End')

ax0.axis(aa)
ax1.axis(aa)
ax2.axis(aa)
ax3.axis(aa)
pfun.dar(ax0)
pfun.dar(ax1)
pfun.dar(ax2)
pfun.dar(ax3)

fn_fig = Ldir['LOo'] / 'plots' / 'tplot_trackmap_mechanism2.png'
plt.savefig(fn_fig)

dsr.close()
dsg.close()