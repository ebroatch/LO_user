"""
Plot results of a particle tracking experiment.
"""
from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()
from lo_user_tools import llxyfun

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
# get track data
fn = in_dir0 / exp_name / rel
dsr = xr.open_dataset(fn, decode_times=False)
NT, NP = dsr.lon.shape
lon_vals = dsr.lon.values
lat_vals = dsr.lat.values
z_start = dsr.z.values[0, :] #starting depth of the particles
lon_start = dsr.lon.values[0,:] #should we add newaxis?
lat_start = dsr.lat.values[0,:]
time_hours = dsr.Time.values
dsr.close()
print('got track data')
# get some grid fields
fng = in_dir0 / exp_name / 'grid.nc'
dsg = xr.open_dataset(fng)
lonp, latp = pfun.get_plon_plat(dsg.lon_rho.values, dsg.lat_rho.values)
hh = dsg.h.values
maskr = dsg.mask_rho.values
dsg.close()
print('got grid data')

#define longitude for ends of sills
#NOTE: THIS IS ONLY FOR THE 20KM SILL CURRENTLY!!!
#add some if statements to change sillland depending on grid name
sillsea = llxyfun.x2lon(40e3,0,45)
sillland = llxyfun.x2lon(60e3,0,45)

#get boolean arrays
lon_start_in = lon_start >= sillland #boolean array for particles starting in the inner basin
# lon_start_insill = lon_start >= sillsea #boolean array for particles starting in the inner basin and sill
lon_in = lon_vals >= sillland #boolean array of particles in the inner basin over time
# lon_insill = lon_vals >= sillsea #boolean array of particles in the inner basin and sill over time
lon_est = lon_vals >= 0 #boolean array of particles in the whole estuary over time
print('got boolean lon arrays\n')

#get strict residence times
tmax = time_hours[-1]
#inner basin
rt_strict_in = np.argmin(lon_in,axis=0).astype('float') #first time the particle is outside the inner basin
rt_strict_in = np.where(rt_strict_in==0, tmax+1, rt_strict_in) #replace 0 with tmax+1, this sets the particles that never leave or ones that were released outside to tmax+1
rt_strict_in = rt_strict_in * lon_start_in #this resets the particles that are not released in the inner basin to zero (necessary?)
rt_strict_in = np.where(rt_strict_in==0, np.nan, rt_strict_in) #this sets particles that are not released in the inner basin or sill to nan
rt_strict_days_in = rt_strict_in/24 #convert to days
# #inner basin + sill
# rt_strict_insill = np.argmin(lon_insill,axis=0).astype('float')
# rt_strict_insill = np.where(rt_strict_insill==0, tmax+1, rt_strict_insill) #replace 0 with tmax+1, this sets the particles that never leave or ones that were released outside to tmax+1
# rt_strict_insill = rt_strict_insill * lon_start_insill #this resets the particles that are not released in the inner basin or sill to zero (necessary?)
# rt_strict_insill = np.where(rt_strict_insill==0, np.nan, rt_strict_insill) #this sets particles that are not released in the inner basin or sill to nan
# rt_strict_days_insill = rt_strict_insill/24 #convert to days
#whole estuary
rt_strict_est = np.argmin(lon_est,axis=0).astype('float')
rt_strict_est = np.where(rt_strict_est==0, tmax+1, rt_strict_est) #replace 0 with tmax+1 (every particle is in the estuary at t=0, these are particles that never leave estuary)
rt_strict_days_est = rt_strict_est/24 #convert to days
print('got residence times\n')


# subsample output for plotting
# npmax = 600 # max number of points to plot
# step = max(1,int(np.floor(NP/npmax)))
# lon = dsr.lon.values[:,::step]
# lat = dsr.lat.values[:,::step]

# select particles originating in a single layer
#depths = np.array([-12.5, -37.5, -62.5, -87.5, -112.5, -137.5, -162.5, -187.5])
# depths = np.array([-12.5, -62.5, -137.5, -187.5])
depths = np.array([-12.5, -112.5, -187.5])

#set up plot
plt.close('all')
#fig, axs = plt.subplots(2,4, sharex=True, sharey=True, figsize=(15,10))
#fig, axs = plt.subplots(2,1, sharex=True, sharey=True)
fig = plt.figure(figsize=(12,12))
gs = fig.add_gridspec(nrows=3,ncols=3, width_ratios=[10,5,1], height_ratios=[1,1,1])
ax1 = fig.add_subplot(gs[0,0]) 
ax2 = fig.add_subplot(gs[1,0])
ax3 = fig.add_subplot(gs[2,0])
ax1b = fig.add_subplot(gs[0,1]) 
ax2b = fig.add_subplot(gs[1,1])
ax3b = fig.add_subplot(gs[2,1])
# ax4 = fig.add_subplot(gs[3,0])
axc = fig.add_subplot(gs[:,2])
# axs=[ax1,ax2,ax3,ax4,ax5]
# axs=[ax1,ax2,ax3,axc]
axs=[ax1,ax2,ax3]
axsb=[ax1b,ax2b,ax3b]


for j in range(len(depths)):
    depth = depths[j]
    ax=axs[j]
    axb=axsb[j]

    #sort the residence times based on their starting layer
    rt_strict_days_in_depth = np.where((z_start>(depth-5)) & (z_start<(depth+5)),rt_strict_days_in,np.nan) #set all particles starting outside the depth layer to nan
    rt_strict_days_est_depth = np.where((z_start>(depth-5)) & (z_start<(depth+5)),rt_strict_days_est,np.nan) #set all particles starting outside the depth layer to nan
    # lon = dsr.lon.where((dsr.z.sel(Time=0)>(depth-5)) & (dsr.z.sel(Time=0)<(depth+5)),drop=True).values
    # lat = dsr.lat.where((dsr.z.sel(Time=0)>(depth-5)) & (dsr.z.sel(Time=0)<(depth+5)),drop=True).values

    #bin the results in 2d so that we can make a map
    # lont0 = lon[0,:] #starting lat and lon which we will use to bin the residence times
    # latt0 = lat[0,:]
    lonbin=lonp[0,:] #use the lonp and latp values as bin edges (lon_rho are centered within them)
    latbin=latp[:,0]
    # ret = stats.binned_statistic_2d(lon_start,lat_start,rtd,'mean',bins=[lonbin,latbin])
    ret_in = stats.binned_statistic_2d(lon_start,lat_start,rt_strict_days_in_depth,np.nanmean,bins=[lonbin,latbin])
    rt_in_map = ret_in.statistic
    ret_est = stats.binned_statistic_2d(lon_start,lat_start,rt_strict_days_est_depth,np.nanmean,bins=[lonbin,latbin])
    rt_est_map = ret_est.statistic

    #plot
    # cs=ax.pcolormesh(lonbin,latbin,np.transpose(rtdmap),vmin=0,vmax=(tmax-1)/24,cmap=cm.matter)
    cs=ax.pcolormesh(lonbin,latbin,np.transpose(rt_est_map),vmin=0,vmax=tmax/24,cmap=cm.matter)
    csb=axb.pcolormesh(lonbin,latbin,np.transpose(rt_in_map),vmin=0,vmax=tmax/24,cmap=cm.matter)
    
    aa = [0,1.15,44.95,45.05]
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(str(depth)+'m depth estuary residence time')

    aab = [sillland,1.15,44.95,45.05]
    axb.axis(aab)
    pfun.dar(axb)
    axb.set_xlabel('Longitude')
    axb.set_ylabel('Latitude')
    axb.set_title(str(depth)+'m depth inner basin residence time')


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

fig.colorbar(cs,cax=axc,extend='max',label='Residence time [days]')
plt.suptitle('Residence times at different depths')

fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtmap_tracker2.png'
plt.savefig(fn_fig)
#plt.show()
#pfun.end_plot()
plt.close()

dsr.close()
dsg.close()