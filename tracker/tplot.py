"""
Plot results of a particle tracking experiment.
"""


import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()

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

# # subsample output for plotting #SKIP SUBSAMPLING
npmax = 600 # max number of points to plot
step = max(1,int(np.floor(NP/npmax)))

lon = dsr.lon.values[:,::step]
lat = dsr.lat.values[:,::step]
# lon = dsr.lon.values #use this for no subsampling
# lat = dsr.lat.values

# make a mask that is False from the time a particle first leaves the domain
# and onwards
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
pfun.start_plot(figsize=(14,8))
fig = plt.figure()

# MAP
# set domain limits
if False:
    # plot full domain
    aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
else:
    # automatically plot region of particles, with padding
    pad = .02
    aa = [np.nanmin(lon) - pad, np.nanmax(lon) + pad,
    np.nanmin(lat) - pad, np.nanmax(lat) + pad]
    
#ax = fig.add_subplot(121)
ax = fig.add_subplot(111)
zm = -np.ma.masked_where(maskr==0, hh)
plt.pcolormesh(lonp, latp, zm, vmin=-200, vmax=0,
    cmap='terrain', alpha=.25)
ax.axis(aa)
pfun.dar(ax)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_title(exp_name.strip('/'))
# add the tracks (packed [time, particle])
# regular spaghetti plots
ax.plot(lon, lat, '-k', linewidth=.2) #DON'T PLOT LINES FOR NOW #comment out to skip lines
# ax.plot(lon[0,:], lat[0,:], 'og', alpha=.3)
# ax.plot(lon[-1,:], lat[-1,:], 'or', alpha=.3)
ax.plot(lon[0,:], lat[0,:], '.g', alpha=.3, markeredgecolor='none')
ax.plot(lon[-1,:], lat[-1,:], '.r', alpha=.3, markeredgecolor='none')

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
fn_fig = Ldir['LOo'] / 'plots' / 'tplot_subsample.png'
plt.savefig(fn_fig)
pfun.end_plot()

# #PLOTTING - HISTOGRAMS
# fig, axs = plt.subplots(5,1,sharex=True)
# for j in range(5):
#     hour=j*180
#     axs[j].set_title('t='+str(hour)+'h')
#     axs[j].hist(dsr['lon'].sel(Time=hour),bins=20,alpha=0.5)
#     #axs[j].set_ylim(0, 30)

# #PLOTTING - HISTOGRAMS
fig, axs = plt.subplots(5,1,sharex=True, figsize=(10,10))
for j in range(5):
    hour=j*720
    axs[j].set_title('t='+str(hour)+'h')
    axs[j].hist(dsr['lon'].sel(Time=hour),bins=20,range=(-0.8,1.2),alpha=0.5)
#     #axs[j].set_ylim(0, 30)
fn_fig2 = Ldir['LOo'] / 'plots' / 'tplot_hist.png'
plt.savefig(fn_fig2)
# plt.show()

dsr.close()
dsg.close()