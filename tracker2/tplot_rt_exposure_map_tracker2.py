"""
Plot results of a particle tracking experiment.
"""
from lo_tools import Lfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()
from lo_user_tools import llxyfun

import matplotlib.pyplot as plt
import matplotlib.colors
import xarray as xr
import numpy as np
from scipy import stats
from cmocean import cm

def dar2(ax):
    """
    Fixes the plot aspect ratio to be locally Cartesian.
    AND anchor to right side
    """
    yl = ax.get_ylim()
    yav = (yl[0] + yl[1])/2
    ax.set_aspect(1/np.cos(np.pi*yav/180),anchor='E')

#Choose sill and set parameters accordingly
sill_choice = input("Enter sill length choice [km]: ")
sill_choice = int(sill_choice)

if sill_choice==5:
    sillsea = llxyfun.x2lon(40e3,0,45)
    sillland = llxyfun.x2lon(45e3,0,45)
    # lonlim = 1.1
    estlenkm = 85
    # estlenlon = llxyfun.x2lon(estlenkm*1e3,0,45)
    lonliminner = estlenlon+0.1
    silllenlabel = '5km'
    fn = '/data1/ebroatch/LO_output/tracks2/sill5km_t0_xa0/sill5kmest_3d/release_2020.09.01.nc'
    fng = '/data1/ebroatch/LO_output/tracks2/sill5km_t0_xa0/sill5kmest_3d/grid.nc'
elif sill_choice==10:
    sillsea = llxyfun.x2lon(40e3,0,45)
    sillland = llxyfun.x2lon(50e3,0,45)
    # lonlim = 1.2
    estlenkm = 90
    estlenlon = llxyfun.x2lon(estlenkm*1e3,0,45)
    lonliminner = estlenlon+0.1
    silllenlabel = '10km'
    fn = '/data1/ebroatch/LO_output/tracks2/sill10km_t2_xa0/sill10kmest_3d/release_2020.09.01.nc'
    fng = '/data1/ebroatch/LO_output/tracks2/sill10km_t2_xa0/sill10kmest_3d/grid.nc'
elif sill_choice==20:
    sillsea = llxyfun.x2lon(40e3,0,45)
    sillland = llxyfun.x2lon(60e3,0,45)
    # lonlim = 1.3
    estlenkm = 100
    estlenlon = llxyfun.x2lon(estlenkm*1e3,0,45)
    lonliminner = estlenlon+0.1
    silllenlabel = '20km'
    fn = '/data1/ebroatch/LO_output/tracks2/sill20kmdeep_t2_xa0/sill20kmdeepest_3d/release_2020.09.01.nc'
    fng = '/data1/ebroatch/LO_output/tracks2/sill20kmdeep_t2_xa0/sill20kmdeepest_3d/grid.nc'
elif sill_choice==40:
    sillsea = llxyfun.x2lon(40e3,0,45)
    sillland = llxyfun.x2lon(80e3,0,45)
    # lonlim = 1.6
    estlenkm = 120
    estlenlon = llxyfun.x2lon(estlenkm*1e3,0,45)
    lonliminner = estlenlon+0.1
    silllenlabel = '40km'
    fn = '/data1/ebroatch/LO_output/tracks2/sill40km_t2_xa0/sill40kmest_3d/release_2020.09.01.nc'
    fng = '/data1/ebroatch/LO_output/tracks2/sill40km_t2_xa0/sill40kmest_3d/grid.nc'
elif sill_choice==80:
    sillsea = llxyfun.x2lon(40e3,0,45)
    sillland = llxyfun.x2lon(120e3,0,45)
    # lonlim = 2.1
    estlenkm = 160
    estlenlon = llxyfun.x2lon(estlenkm*1e3,0,45)
    lonliminner = estlenlon+0.1
    silllenlabel = '80km'
    fn = '/data1/ebroatch/LO_output/tracks2/sill80km_t2_xa0/sill80kmest_3d/release_2020.09.01.nc'
    fng = '/data1/ebroatch/LO_output/tracks2/sill80km_t2_xa0/sill80kmest_3d/grid.nc'
else:
    print('Sill length must be 5, 10, 20, 40, or 80')

# # Choose an experiment and release to plot.
# in_dir0 = Ldir['LOo'] / 'tracks'
# exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
#     itext='** Choose experiment from list **', last=False)
# rel = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
#     itext='** Choose item from list **', last=False)

# get Datasets
# get track data
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
dsg = xr.open_dataset(fng)
lonp, latp = pfun.get_plon_plat(dsg.lon_rho.values, dsg.lat_rho.values)
hh = dsg.h.values
maskr = dsg.mask_rho.values
dsg.close()
print('got grid data')

#get boolean arrays
lon_est = lon_vals >= 0 #boolean array of particles in the whole estuary over time
lon_in = lon_vals >= sillland #boolean array of particles in the inner basin over time
lon_start_in = lon_start >= sillland #boolean array for particles starting in the inner basin
# lon_start_insill = lon_start >= sillsea #boolean array for particles starting in the inner basin and sill
# lon_insill = lon_vals >= sillsea #boolean array of particles in the inner basin and sill over time
print('got boolean lon arrays\n')

#get strict residence times
tmax = time_hours[-1]
#whole estuary
rt_strict_est = np.argmin(lon_est,axis=0).astype('float')
rt_strict_est = np.where(rt_strict_est==0, tmax+1, rt_strict_est) #replace 0 with tmax+1 (every particle is in the estuary at t=0, these are particles that never leave estuary)
rt_strict_days_est = rt_strict_est/24 #convert to days
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
print('got residence times\n')

print('average estuary residence time [days]')
print(np.nanmean(rt_strict_days_est))
print('% of particles with residence time >120d')
print(100*np.count_nonzero(rt_strict_est==tmax+1)/np.count_nonzero(~np.isnan(rt_strict_est)))
print('average inner basin residence time [days]')
print(np.nanmean(rt_strict_days_in))
print('% of particles with residence time >120d')
print(100*np.count_nonzero(rt_strict_in==tmax+1)/np.count_nonzero(~np.isnan(rt_strict_in)))


#get exposure times
exposuret_est = np.sum(lon_est,axis=0)
exposuret_days_est = exposuret_est/24
exposuret_in = np.sum(lon_in,axis=0)
exposuret_in = exposuret_in * lon_start_in #this resets the particles that are not released in the inner basin to zero (necessary?)
exposuret_in = np.where(exposuret_in==0, np.nan, exposuret_in) #this sets particles that are not released in the inner basin or sill to nan
exposuret_days_in = exposuret_in/24

print('average estuary exposure time [days]')
print(np.nanmean(exposuret_days_est))
print('average inner basin exposure time [days]')
print(np.nanmean(exposuret_days_in))

# subsample output for plotting
# npmax = 600 # max number of points to plot
# step = max(1,int(np.floor(NP/npmax)))
# lon = dsr.lon.values[:,::step]
# lat = dsr.lat.values[:,::step]

# select particles originating in a single layer
#depths = np.array([-12.5, -37.5, -62.5, -87.5, -112.5, -137.5, -162.5, -187.5])
# depths = np.array([-12.5, -62.5, -137.5, -187.5])
# depths = np.array([-12.5, -112.5, -187.5])
depths = np.array([-12.5, -187.5]) #top and bottom layers only
depthstr = ['12.5m','187.5m']

#set up plot
plt.close('all')
#fig, axs = plt.subplots(2,4, sharex=True, sharey=True, figsize=(15,10))
#fig, axs = plt.subplots(2,1, sharex=True, sharey=True)
# wr1= 8*estlenkm/40
# wr1=8*((lonlim)/(lonlim-sillland))
lonlim = 2.15
wr1=8*((lonlim)/(lonliminner-sillland))
fig = plt.figure(figsize=(20,7))
# fig = plt.figure(figsize=(12,6))
gs = fig.add_gridspec(nrows=4,ncols=3, width_ratios=[wr1,8,1], height_ratios=[1,1,1,1])
ax1 = fig.add_subplot(gs[0,0])
ax1b = fig.add_subplot(gs[1,0])  
ax2 = fig.add_subplot(gs[2,0])
ax2b = fig.add_subplot(gs[3,0]) 
ax1c = fig.add_subplot(gs[0,1])
ax1d = fig.add_subplot(gs[1,1]) 
ax2c = fig.add_subplot(gs[2,1])
ax2d = fig.add_subplot(gs[3,1])
axcb = fig.add_subplot(gs[:,2])
# axs=[ax1,ax2,ax3,ax4,ax5]
# axs=[ax1,ax2,ax3,axc]
axs=[ax1,ax2]
axsb=[ax1b,ax2b]
axsc=[ax1c,ax2c]
axsd=[ax1d,ax2d]


for j in range(len(depths)):
    depth = depths[j]
    ax=axs[j]
    axb=axsb[j]
    axc=axsc[j]
    axd=axsd[j]

    #sort the residence times based on their starting layer
    rt_strict_days_est_depth = np.where((z_start>(depth-5)) & (z_start<(depth+5)),rt_strict_days_est,np.nan) #set all particles starting outside the depth layer to nan
    rt_strict_days_in_depth = np.where((z_start>(depth-5)) & (z_start<(depth+5)),rt_strict_days_in,np.nan) #set all particles starting outside the depth layer to nan
    exposuret_days_est_depth = np.where((z_start>(depth-5)) & (z_start<(depth+5)),exposuret_days_est,np.nan) #set all particles starting outside the depth layer to nan
    exposuret_days_in_depth = np.where((z_start>(depth-5)) & (z_start<(depth+5)),exposuret_days_in,np.nan) #set all particles starting outside the depth layer to nan
    # lon = dsr.lon.where((dsr.z.sel(Time=0)>(depth-5)) & (dsr.z.sel(Time=0)<(depth+5)),drop=True).values
    # lat = dsr.lat.where((dsr.z.sel(Time=0)>(depth-5)) & (dsr.z.sel(Time=0)<(depth+5)),drop=True).values


    #bin the results in 2d so that we can make a map
    # lont0 = lon[0,:] #starting lat and lon which we will use to bin the residence times
    # latt0 = lat[0,:]
    lonbin=lonp[0,:] #use the lonp and latp values as bin edges (lon_rho are centered within them)
    latbin=latp[:,0]
    # ret = stats.binned_statistic_2d(lon_start,lat_start,rtd,'mean',bins=[lonbin,latbin])
    ret_est = stats.binned_statistic_2d(lon_start,lat_start,rt_strict_days_est_depth,np.nanmean,bins=[lonbin,latbin])
    rt_est_map = ret_est.statistic
    ret_in = stats.binned_statistic_2d(lon_start,lat_start,rt_strict_days_in_depth,np.nanmean,bins=[lonbin,latbin])
    rt_in_map = ret_in.statistic
    ret_est_exposure = stats.binned_statistic_2d(lon_start,lat_start,exposuret_days_est_depth,np.nanmean,bins=[lonbin,latbin])
    exposuret_est_map = ret_est_exposure.statistic
    ret_in_exposure = stats.binned_statistic_2d(lon_start,lat_start,exposuret_days_in_depth,np.nanmean,bins=[lonbin,latbin])
    exposuret_in_map = ret_in_exposure.statistic


    #plot
    levels = [0,10,20,30,40,50,60,70,80,90,100,110,120]
    cmap = plt.colormaps['turbo']
    norm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N, extend='max')
    # cs=ax.pcolormesh(lonbin,latbin,np.transpose(rtdmap),vmin=0,vmax=(tmax-1)/24,cmap=cm.matter)
    cs=ax.pcolormesh(lonbin,latbin,np.transpose(rt_est_map),cmap='turbo',norm=norm)
    csb=axb.pcolormesh(lonbin,latbin,np.transpose(exposuret_est_map),cmap='turbo',norm=norm)
    csc=axc.pcolormesh(lonbin,latbin,np.transpose(rt_in_map),cmap='turbo',norm=norm)
    csd=axd.pcolormesh(lonbin,latbin,np.transpose(exposuret_in_map),cmap='turbo',norm=norm)
    # cs=ax.contourf((lonbin[1:]+lonbin[:-1])/2,(latbin[1:]+latbin[:-1])/2,np.transpose(rt_est_map),cmap='turbo',levels=[0,10,20,30,40,50,60,70,80,90,100,110,120],extend='max') #try with contour
    # csb=axb.contourf((lonbin[1:]+lonbin[:-1])/2,(latbin[1:]+latbin[:-1])/2,np.transpose(rt_in_map),cmap='turbo',levels=[0,10,20,30,40,50,60,70,80,90,100,110,120],extend='max')
    
    #set limits
    aa = [0,lonlim,44.95,45.05]
    ax.axis(aa)
    pfun.dar(ax)
    # ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    axb.axis(aa)
    pfun.dar(axb)
    axb.set_ylabel('Latitude')
    if j==(len(depths)-1):
        axb.set_xlabel('Longitude')
    ax.set_title('Whole estuary residence time ('+depthstr[j]+' depth)')
    axb.set_title('Whole estuary exposure time ('+depthstr[j]+' depth)')

    aac = [sillland,lonliminner,44.95,45.05]
    axc.axis(aac)
    pfun.dar(axc)
    # axc.set_xlabel('Longitude')
    # axc.set_ylabel('Latitude')
    axd.axis(aac)
    pfun.dar(axd)
    # axd.set_ylabel('Latitude')   
    if j==(len(depths)-1):
        axd.set_xlabel('Longitude')
    axc.set_title('Inner basin residence time ('+depthstr[j]+' depth)')
    axd.set_title('Inner basin exposure time ('+depthstr[j]+' depth)')
    axc.yaxis.set_tick_params(labelleft=True)
    axd.yaxis.set_tick_params(labelleft=True)


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

fig.colorbar(cs,cax=axcb,extend='max',label='Time scale [days]')
# plt.suptitle('Residence and exposure times at different depths for '+silllenlabel+' sill model')
#ADD SILL LENGTH!!!!
axs[0].text(0, 45.051, silllenlabel+' sill model', horizontalalignment='left', verticalalignment='bottom')

# axs[0].text(estlenlon+0.05, 45, 'A', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
# axs[1].text(estlenlon+0.05, 45, 'B', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
# axsb[0].text(estlenlon+0.05, 45, 'C', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
# axsb[1].text(estlenlon+0.05, 45, 'D', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
# axsc[0].text(estlenlon+0.05, 45, 'E', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
# axsc[1].text(estlenlon+0.05, 45, 'F', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
# axsd[0].text(estlenlon+0.05, 45, 'G', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
# axsd[1].text(estlenlon+0.05, 45, 'H', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
axs[0].text(2.1, 45, 'A', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
axsc[0].text(estlenlon+0.05, 45, 'B', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
axsb[0].text(2.1, 45, 'C', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
axsd[0].text(estlenlon+0.05, 45, 'D', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
axs[1].text(2.1, 45, 'E', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
axsc[1].text(estlenlon+0.05, 45, 'F', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
axsb[1].text(2.1, 45, 'G', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')
axsd[1].text(estlenlon+0.05, 45, 'H', horizontalalignment='center', verticalalignment='center', fontsize=14, fontweight='bold')

fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rt_exposure_map_tracker2.png'
plt.savefig(fn_fig)
#plt.show()
#pfun.end_plot()
plt.close()

dsr.close()
dsg.close()