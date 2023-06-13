"""
Plot results of a particle tracking experiment.
"""
from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from scipy import stats
from cmocean import cm

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
#depths = np.array([-12.5, -37.5, -62.5, -87.5, -112.5, -137.5, -162.5, -187.5])
depths = np.array([-12.5, -62.5, -137.5, -187.5])
plt.close('all')
#fig, axs = plt.subplots(2,4, sharex=True, sharey=True, figsize=(15,10))
#fig, axs = plt.subplots(2,1, sharex=True, sharey=True)
fig = plt.figure(figsize=(8,12))
gs = fig.add_gridspec(nrows=4,ncols=2, width_ratios=[10,1], height_ratios=[1,1,1,1])
ax1 = fig.add_subplot(gs[0,0]) 
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[2,0])
ax4 = fig.add_subplot(gs[3,0])
ax5 = fig.add_subplot(gs[:,1])
axs=[ax1,ax2,ax3,ax4,ax5]


for j in range(4):
    depth = depths[j]
    ax=axs[j]
    lon = dsr.lon.where((dsr.z.sel(Time=0)>(depth-5)) & (dsr.z.sel(Time=0)<(depth+5)),drop=True).values
    lat = dsr.lat.where((dsr.z.sel(Time=0)>(depth-5)) & (dsr.z.sel(Time=0)<(depth+5)),drop=True).values
    lont0 = lon[0,:]
    latt0 = lat[0,:]
    tmax = lon.shape[0]
    rt = np.argmax(lon<0,axis=0).astype('float')
    rt = np.where(rt==0, tmax, rt) #replace 0 with tmax (these particles never leave estuary - it's really just greater than tmax)
    rtd = rt/24 #convert to days
    lonbin=lonp[0,:]
    latbin=latp[:,0]
    ret = stats.binned_statistic_2d(lont0,latt0,rtd,'mean',bins=[lonbin,latbin])
    rtdmap = ret.statistic
    cs=ax.pcolormesh(lonbin,latbin,np.transpose(rtdmap),vmin=0,vmax=(tmax-1)/24,cmap=cm.matter)
    
    aa = [0,1.15,44.95,45.05]
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(str(depth)+'m')


    # # make a mask that is False from the time a particle first leaves the domain
    # # and onwards
    # AA = [dsg.lon_rho.values[0,0], dsg.lon_rho.values[0,-1],
    #         dsg.lat_rho.values[0,0], dsg.lat_rho.values[-1,0]]
    # ib_mask = np.ones(lon.shape, dtype=bool)
    # ib_mask[lon < AA[0]] = False
    # ib_mask[lon > AA[1]] = False
    # ib_mask[lat < AA[2]] = False
    # ib_mask[lat > AA[3]] = False
    # NTS, NPS = lon.shape
    # for pp in range(NPS):
    #     tt = np.argwhere(ib_mask[:,pp]==False)
    #     if len(tt) > 0:
    #         ib_mask[tt[0][0]:, pp] = False

    # # and apply the mask to lon and lat
    # lon[~ib_mask] = np.nan
    # lat[~ib_mask] = np.nan

    # PLOTTING - SPAGHETTI PLOT

    # # MAP
    # # set domain limits
    # if True:
    #     # plot full domain
    #     #aa = [lonp.min(), lonp.max(), latp.min(), latp.max()]
    #     aa = [-0.2,1.2,43.5,46.5]
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
    # ax.set_xlabel('Longitude')
    # ax.set_ylabel('Latitude')
    # ax.set_title(str(depth)+'m')
    # add the tracks (packed [time, particle])
    # regular spaghetti plots
    # ax.plot(lon, lat, '-k', linewidth=.2)
    # ax.plot(lon[0,:], lat[0,:], 'og', alpha=.3)
    # ax.plot(lon[-1,:], lat[-1,:], 'or', alpha=.3)

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
    # ax.plot(lon[0,:], lat[0,:], '.g', alpha=.3, markeredgecolor='none')
    #ax.plot(lon[-1,:], lat[-1,:], '.r', alpha=.3, markeredgecolor='none')

fig.colorbar(cs,cax=axs[-1],extend='max',label='Residence time [days]')
plt.suptitle('Residence times at different depths')

fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtmap.png'
plt.savefig(fn_fig)
#plt.show()
#pfun.end_plot()
plt.close()

dsr.close()
dsg.close()