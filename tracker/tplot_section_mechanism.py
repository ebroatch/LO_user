"""
Plot results of a particle tracking experiment.
"""


import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
from lo_user_tools import llxyfun
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

# subsample output for plotting
# npmax = 600 # max number of points to plot
# step = max(1,int(np.floor(NP/npmax)))
# lon = dsr.lon.values[:,::step]
# lat = dsr.lat.values[:,::step]

# select particles originating in a single layer
plt.close('all')
#fig, ax = plt.subplots(1,1)
fig, [ax,ax2] = plt.subplots(2,1,figsize=(12,10))
seclat = 45
#depths = np.array([-12.5, -37.5, -62.5, -87.5, -112.5, -137.5, -162.5, -187.5])
#depth = -12.5
#sillmid = llxyfun.x2lon(50e3,0,45) #cut out particles starting on sill, use sillsea and sillland instead
sillsea = llxyfun.x2lon(40e3,0,45)
sillland = llxyfun.x2lon(60e3,0,45)
# lon = dsr.lon.where((dsr.lat.sel(Time=0)==seclat) & (dsr.z.sel(Time=0)>(depth-5)) & (dsr.z.sel(Time=0)<(depth+5)) & (dsr.lon.sel(Time=0)<sillmid),drop=True).values
# z = dsr.z.where((dsr.lat.sel(Time=0)==seclat) & (dsr.z.sel(Time=0)>(depth-5)) & (dsr.z.sel(Time=0)<(depth+5)) & (dsr.lon.sel(Time=0)<sillmid),drop=True).values
# lon2 = dsr.lon.where((dsr.lat.sel(Time=0)==seclat) & (dsr.z.sel(Time=0)>(depth-5)) & (dsr.z.sel(Time=0)<(depth+5)) & (dsr.lon.sel(Time=0)>sillmid),drop=True).values
# z2 = dsr.z.where((dsr.lat.sel(Time=0)==seclat) & (dsr.z.sel(Time=0)>(depth-5)) & (dsr.z.sel(Time=0)<(depth+5)) & (dsr.lon.sel(Time=0)>sillmid),drop=True).values

# lon = dsr.lon.where((dsr.lat.sel(Time=0)==seclat) & (dsr.lon.sel(Time=0)<sillsea),drop=True).values
# z = dsr.z.where((dsr.lat.sel(Time=0)==seclat) & (dsr.lon.sel(Time=0)<sillsea),drop=True).values
# lon2 = dsr.lon.where((dsr.lat.sel(Time=0)==seclat) & (dsr.lon.sel(Time=0)>sillland),drop=True).values
# z2 = dsr.z.where((dsr.lat.sel(Time=0)==seclat) & (dsr.lon.sel(Time=0)>sillland),drop=True).values

lon = dsr.lon.where((dsr.lon.sel(Time=0)<sillland) & (dsr.lon.sel(Time=0)>sillsea),drop=True).values
z = dsr.z.where((dsr.lon.sel(Time=0)<sillland) & (dsr.lon.sel(Time=0)>sillsea),drop=True).values
lon2 = dsr.lon.where((dsr.lon.sel(Time=12)<sillland) & (dsr.lon.sel(Time=12)>sillsea),drop=True).values
z2 = dsr.z.where((dsr.lon.sel(Time=12)<sillland) & (dsr.lon.sel(Time=12)>sillsea),drop=True).values


hg45 = dsg.h.where(dsg.lat_rho==seclat, drop=True).values.squeeze() #profile of bottom bathymetry
long45 = dsg.lon_rho.where(dsg.lat_rho==seclat, drop=True).values.squeeze()

# make a mask that is False from the time a particle first leaves the domain
# and onwards
AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
        dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
ib_mask = np.ones(lon.shape, dtype=bool)
ib_mask[lon < AA[0]] = False
ib_mask[lon > AA[1]] = False
# ib_mask[lat < AA[2]] = False
# ib_mask[lat > AA[3]] = False
NTS, NPS = lon.shape
for pp in range(NPS):
    tt = np.argwhere(ib_mask[:,pp]==False)
    if len(tt) > 0:
        ib_mask[tt[0][0]:, pp] = False

ib_mask2 = np.ones(lon2.shape, dtype=bool)
ib_mask2[lon2 < AA[0]] = False
ib_mask2[lon2 > AA[1]] = False
# ib_mask[lat < AA[2]] = False
# ib_mask[lat > AA[3]] = False
NTS2, NPS2 = lon2.shape
for pp in range(NPS2):
    tt2 = np.argwhere(ib_mask2[:,pp]==False)
    if len(tt2) > 0:
        ib_mask2[tt2[0][0]:, pp] = False

# and apply the mask to lon and lat
lon[~ib_mask] = np.nan
# lat[~ib_mask] = np.nan
lon2[~ib_mask2] = np.nan

# PLOTTING - SPAGHETTI PLOT

# # MAP
# # set domain limits
# if True:
#     # plot full domain
#     #aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
#     aa = [-2,1.5,43.5,46.5]
# else:
#     # automatically plot region of particles, with padding
#     pad = .02
#     aa = [np.nanmin(lon) - pad, np.nanmax(lon) + pad,
#     np.nanmin(lat) - pad, np.nanmax(lat) + pad]
    
# zm = -np.ma.masked_where(maskr==0, hh)
# ax.pcolormesh(lonp, latp, zm, vmin=-200, vmax=0,
#     cmap='terrain', alpha=.25)
# ax.axis(aa)
# pfun.dar(ax)
ax2.set_xlabel('Longitude')
ax.set_ylabel('z')
ax2.set_ylabel('z')
ax.set_title('12h tracks starting on sill')
ax2.set_title('12h tracks ending on sill')

# add the tracks (packed [time, particle])
# regular spaghetti plots
#step = 40 #step for subsampling lines
step = 1 #no subsmaple
ax.plot(lon[:13,::step], z[:13,::step], '-', color='tab:cyan', linewidth=.1, label='Track') #plot hours 0 to 12
ax.plot(lon[12,:], z[12,:], '.', color='tab:blue', label='End')
# ax.plot(lon[0,:], z[0,:], '.', color='tab:blue', label='Start')
# ax.plot(lon[12,:], z[12,:], '*', color='b', label='End') #end point is after 12 hours
ax2.plot(lon2[:13,::step], z2[:13,::step], '-', color='tab:pink', linewidth=.1, label='Track') #plot hours 0 to 12
ax2.plot(lon2[0,:], z2[0,:], '.', color='m', label='Start')
#ax2.plot(lon2[12,:], z2[12,:], '*', color='r', label='End')
# ax.plot(lon[0,:], z[0,:], '.g', alpha=.3, markeredgecolor='none')
# ax.plot(lon[-1,:], z[-1,:], '.r', alpha=.3, markeredgecolor='none')
# ax.legend()
# ax2.legend()

# add the bottom bathymetry
ax.plot(long45,-hg45,'-k')
ax2.plot(long45,-hg45,'-k')
ax.plot([0,0],[-200,0],'--k')
ax2.plot([0,0],[-200,0],'--k')
ax.plot([-0.2,0],[-15.72,0],':k')
ax2.plot([-0.2,0],[-15.72,0],':k')

# axis limits
ax.set_ylim(-205,5)
#ax.set_xlim(-0.2,1.3)
ax.set_xlim(0,1.3)
ax2.set_ylim(-205,5)
#ax2.set_xlim(-0.2,1.3)
ax2.set_xlim(0,1.3)
# ax.legend()
# ax2.legend()

# time series
# td = (ot_vec - ot_vec[0])/86400
# tv_list = ['z', 'u', 'v']
# #tv_list = ['u', 'v', 'lon', 'lat']
# ntv = len(tv_list)
# for ii in range(ntv):
#     tv = tv_list[ii]
#     NC = 2
#     ax = fig.add_subplot(ntv,NC, (ii+1)*NC)
#     v = dsr[tv].values[:,::step]
#     v[~ib_mask] = np.nan
#     ax.plot(td, v, lw=.5, alpha=.2)
#     ax.text(.05, .05, tv, fontweight='bold', transform=ax.transAxes)
#     if ii == ntv-1:
#         ax.set_xlabel('Time (days)')

#plt.suptitle('exp_name.strip('/')')

#plt.show()
pfun.end_plot()

#plt.show()
fn_fig = Ldir['LOo'] / 'plots' / 'tplot_section_mechanism.png'
plt.savefig(fn_fig)
plt.close()
dsr.close()
dsg.close()