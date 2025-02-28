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
in_dir0 = Ldir['LOo'] / 'tracks2'
gtx_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    itext='** Choose run from list **', last=False)
exp_name = Lfun.choose_item(in_dir0 / gtx_name, tag='', exclude_tag='.csv',
    itext='** Choose experiment from list **', last=False)
rel = Lfun.choose_item(in_dir0 / gtx_name / exp_name, tag='.nc', exclude_tag='grid',
    itext='** Choose item from list **', last=False)

# get Datasets
fn = in_dir0 / gtx_name / exp_name / rel
fng = in_dir0 / gtx_name / exp_name / 'grid.nc'
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

# # subsample output for plotting #SKIP SUBSAMPLING
npmax = 150 # max number of points to plot
step = max(1,int(np.floor(NP/npmax))) #use same for both assuming approx equal number of particles in each basin
# sillmid = llxyfun.x2lon(44e3,0,45)
# lon1 = dsr.lon.where((dsr.lon.sel(Time=0)<sillmid),drop=True).values[:,::step]
# lon2 = dsr.lon.where((dsr.lon.sel(Time=0)>sillmid),drop=True).values[:,::step]
# lat1 = dsr.lat.where((dsr.lon.sel(Time=0)<sillmid),drop=True).values[:,::step]
# lat2 = dsr.lat.where((dsr.lon.sel(Time=0)>sillmid),drop=True).values[:,::step]
sillsea = llxyfun.x2lon(40e3,0,45)
# sillland = llxyfun.x2lon(60e3,0,45)
if gtx_name.split('_')[0]=='sill5km':
    sillland = llxyfun.x2lon(45e3,0,45)
    xlonlim=1.1
elif gtx_name.split('_')[0]=='sill10km':
    sillland = llxyfun.x2lon(50e3,0,45)
    xlonlim=1.2
#elif grid=='20kmdeep':
elif gtx_name.split('_')[0]=='sill20kmdeep':
    sillland = llxyfun.x2lon(60e3,0,45)
    xlonlim=1.3
elif gtx_name.split('_')[0]=='sill40km':
    sillland = llxyfun.x2lon(80e3,0,45)
    xlonlim=1.6
# elif grid=='80km':
elif gtx_name.split('_')[0]=='sill80km':
    sillland = llxyfun.x2lon(120e3,0,45)
    xlonlim=2.1

# lon1 = dsr.lon.where((dsr.lon.sel(Time=0)<sillsea),drop=True).values[:,::step]
# print('got lon for outer basin')
# lon2 = dsr.lon.where((dsr.lon.sel(Time=0)>sillland),drop=True).values[:,::step]
# print('got lon for inner basin')
# lat1 = dsr.lat.where((dsr.lon.sel(Time=0)<sillsea),drop=True).values[:,::step]
# print('got lat for outer basin')
# lat2 = dsr.lat.where((dsr.lon.sel(Time=0)>sillland),drop=True).values[:,::step]
# print('got lat for inner basin')


lon = dsr.lon.values[:,::step]
lat = dsr.lat.values[:,::step]
# lon = dsr.lon.values #use this for no subsampling
# lat = dsr.lat.values

# # make a mask that is False from the time a particle first leaves the domain
# # and onwards (outer basin particles)
# AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
#         dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
# ib_mask = np.ones(lon1.shape, dtype=bool)
# ib_mask[lon1 < AA[0]] = False
# ib_mask[lon1 > AA[1]] = False
# ib_mask[lat1 < AA[2]] = False
# ib_mask[lat1 > AA[3]] = False
# NTS, NPS = lon1.shape
# for pp in range(NPS):
#     tt = np.argwhere(ib_mask[:,pp]==False)
#     if len(tt) > 0:
#         ib_mask[tt[0][0]:, pp] = False

# # and apply the mask to lon and lat
# lon1[~ib_mask] = np.nan
# lat1[~ib_mask] = np.nan

# # make a mask that is False from the time a particle first leaves the domain 
# # and onwards (inner basin particles)
# AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
#         dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
# ib_mask = np.ones(lon2.shape, dtype=bool)
# ib_mask[lon2 < AA[0]] = False
# ib_mask[lon2 > AA[1]] = False
# ib_mask[lat2 < AA[2]] = False
# ib_mask[lat2 > AA[3]] = False
# NTS, NPS = lon2.shape
# for pp in range(NPS):
#     tt = np.argwhere(ib_mask[:,pp]==False)
#     if len(tt) > 0:
#         ib_mask[tt[0][0]:, pp] = False

# # and apply the mask to lon and lat
# lon2[~ib_mask] = np.nan
# lat2[~ib_mask] = np.nan

# make a mask that is False from the time a particle first leaves the domain
# and onwards (outer basin particles)
AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
        dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
ib_mask = np.ones(lon.shape, dtype=bool)
ib_mask[lon < AA[0]] = False
ib_mask[lon > AA[1]] = False
ib_mask[lat < AA[2]] = False
ib_mask[lat > AA[3]] = False
NTS, NPS = lon.shape
for pp in range(NPS):
    tt = np.argwhere(ib_mask[:,pp]==False)
    if len(tt) > 0:
        ib_mask[tt[0][0]:, pp] = False
# and apply the mask to lon and lat
lon[~ib_mask] = np.nan
lat[~ib_mask] = np.nan

# PLOTTING - SPAGHETTI PLOT
plt.close('all')
# pfun.start_plot(figsize=(8,10))
# fig = plt.figure()
fig, [ax,ax2] = plt.subplots(2,1,figsize=(8,10),gridspec_kw={'height_ratios': [8.33,1]})
fig, ax = plt.subplots(1,1,figsize=(10,10))

# MAP
# set domain limits
if False:
    # plot full domain
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
else:
    # automatically plot region of particles, with padding #use outer particles to pick limits
    pad = .02
    aa = [np.nanmin(lon) - pad, np.nanmax(lon) + pad,
    np.nanmin(lat) - pad, np.nanmax(lat) + pad]
    
#ax = fig.add_subplot(121)
# ax = fig.add_subplot(111)
zm = -np.ma.masked_where(maskr==0, hh)
# ax.pcolormesh(lonp, latp, zm, vmin=-200, vmax=20,
#     cmap='bone')
# ax2.pcolormesh(lonp, latp, zm, vmin=-200, vmax=20,
#     cmap='bone')
# plt.pcolormesh(lonp, latp, zm, vmin=-300, vmax=0,
    # cmap='Spectral_r',alpha=0.3)
ax.pcolormesh(lonp, latp, zm, vmin=-300, vmax=0, cmap='Spectral_r')
# ax2.pcolormesh(lonp, latp, zm, vmin=-300, vmax=0, cmap='Spectral_r')
# ax.axis(aa)
ax.set_ylim(43,47)
ax.set_xlim(-4,1.5)
# ax2.set_ylim(44.9,45.1)
# ax2.set_xlim(-0.5,1.5)
pfun.dar(ax)
# pfun.dar(ax2)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
# ax2.set_xlabel('Longitude')
# ax2.set_ylabel('Latitude')
# ax.set_title(exp_name.strip('/'))
# ax.text(.02, .98, exp_name.strip('/'), horizontalalignment='left', verticalalignment='top', transform=ax.transAxes)
ax.text(.02, .02, '20km sill model', horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
# add the tracks (packed [time, particle])
# regular spaghetti plots
# ax.plot(lon1[::2], lat1[::2], '-', color='tab:cyan', linewidth=.2) #DON'T PLOT LINES FOR NOW #comment out to skip lines
# ax.plot(lon2[::2], lat2[::2], '-', color='tab:pink', linewidth=.2)
# ax2.plot(lon1[::2], lat1[::2], '-', color='tab:cyan', linewidth=.2) #DON'T PLOT LINES FOR NOW #comment out to skip lines
# ax2.plot(lon2[::2], lat2[::2], '-', color='tab:pink', linewidth=.2)
ax.plot(lon, lat, '-', color='k', linewidth=.2) #black lines
# ax2.plot(lon, lat, '-', color='k', linewidth=.2)
# ax.plot(lon[0,:], lat[0,:], 'og', alpha=.3)
# ax.plot(lon[-1,:], lat[-1,:], 'or', alpha=.3)
#ax.plot(lon[0,:], lat[0,:], '.g', alpha=.3, markeredgecolor='none')
#ax.plot(lon[-1,:], lat[-1,:], '.r', alpha=.3, markeredgecolor='none')
ax.plot(lon[-1,:], lat[-1,:], '*', color='tab:grey', alpha=1, markeredgecolor='none')
# ax2.plot(lon[-1,:], lat[-1,:], '*', color='tab:grey', alpha=1, markeredgecolor='none')

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
figname = 'tplot_trackmap_tracker2_'+ gtx_name.split('_')[0]+'.png'
fn_fig = Ldir['LOo'] / 'plots' / figname
plt.savefig(fn_fig)
#plt.show()
pfun.end_plot()

dsr.close()
dsg.close()