"""
Plot results of a particle tracking experiment.
"""
from lo_tools import Lfun
from lo_tools import zfun
from lo_tools import plotting_functions as pfun
from lo_user_tools import llxyfun
Ldir = Lfun.Lstart()
import sys

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats

plt.close('all')
# fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(20,6))
# fig, axs = plt.subplots(2,2,figsize=(15,15))
# fig, ax = plt.subplots(1,1,figsize=(15,8))
# fig, axs = plt.subplots(1,3,figsize=(20,8),gridspec_kw={'width_ratios': [4,3,1]})
# fig, axs = plt.subplots(1,3,figsize=(20,8))
fig, axs = plt.subplots(1,2,figsize=(12,6))

for i in range(5):
    # Choose an experiment and release to plot.
    # in_dir0 = Ldir['LOo'] / 'tracks'
    # exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    #     itext='** Choose experiment from list **', last=False)
    # rel = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    #     itext='** Choose item from list **', last=False)


    #for now just use 5,20,80 - can add 10 and 40 if it runs fast
    if i==0:
        silllen = 5
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(45e3,0,45)
        linecolor = 'tab:red'
        linecolor2 = plt.cm.tab20(7)
        silllenlabel = '5km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill5km_t0_xa0/sill5kmest_3d/release_2020.09.01.nc'
        fng = '/data1/ebroatch/LO_output/tracks2/sill5km_t0_xa0/sill5kmest_3d/grid.nc'
    elif i==1:
        silllen = 10
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(50e3,0,45)
        linecolor = 'tab:orange'
        linecolor2 = plt.cm.tab20(3)
        silllenlabel = '10km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill10km_t2_xa0/sill10kmest_3d/release_2020.09.01.nc'
        fng = '/data1/ebroatch/LO_output/tracks2/sill10km_t2_xa0/sill10kmest_3d/grid.nc'
    elif i==2:
        silllen = 20
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(60e3,0,45)
        linecolor = 'tab:green'
        linecolor2 = plt.cm.tab20(5)
        silllenlabel = '20km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill20kmdeep_t2_xa0/sill20kmdeepest_3d/release_2020.09.01.nc'
        fng = '/data1/ebroatch/LO_output/tracks2/sill20kmdeep_t2_xa0/sill20kmdeepest_3d/grid.nc'
    elif i==3:
        silllen = 40
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(80e3,0,45)
        linecolor = 'tab:blue'
        linecolor2 = plt.cm.tab20(1)
        silllenlabel = '40km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill40km_t2_xa0/sill40kmest_3d/release_2020.09.01.nc'
        fng = '/data1/ebroatch/LO_output/tracks2/sill40km_t2_xa0/sill40kmest_3d/grid.nc'
    elif i==4:
        silllen = 80
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(120e3,0,45)
        linecolor = 'tab:purple'
        linecolor2 = plt.cm.tab20(9)
        silllenlabel = '80km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill80km_t2_xa0/sill80kmest_3d/release_2020.09.01.nc'
        fng = '/data1/ebroatch/LO_output/tracks2/sill80km_t2_xa0/sill80kmest_3d/grid.nc'
    
    print(silllenlabel+'\n')

    # get Datasets
    #fn = in_dir0 / exp_name / rel
    #fng = in_dir0 / exp_name / 'grid.nc'
    dsr = xr.open_dataset(fn, decode_times=False)
    dsg = xr.open_dataset(fng)
    lonp, latp = pfun.get_plon_plat(dsg.lon_rho.values, dsg.lat_rho.values) #will use these for the spatial bins
    dsg.close()

    NT, NP = dsr.lon.shape

    # get a list of datetimes
    # ot_vec = dsr.ot.values
    # dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

    # # subsample output for plotting #SKIP SUBSAMPLING
    # npmax = 600 # max number of points to plot
    # step = max(1,int(np.floor(NP/npmax)))

    # lon = dsr.lon.values[:,::step]
    # lat = dsr.lat.values[:,::step]
    # lon = dsr.lon


    lon_vals = dsr.lon.values
    #s_start = dsr.salt.values[0,:]
    z_start = dsr.z.values[0, :] #starting depth of the particles
    lon_start = dsr.lon.values[0,:] #should we add newaxis?
    time_hours = dsr.Time.values
    dsr.close()
    print('got lon_vals and time\n')

    # z_start_lower = z_start < -50 #boolean array for particles starting below sill depth
    # z_start_upper = z_start >= -50 #boolean array for particles starting above sill depth

    # lon_in = lon_vals >= sillland #these are particles in the inner basin at any given point in time #NEED TO CHANGE STRATEGY FOR THESE SMALLER CONTROL VOLUMES, NOT ALL PARTICLES SHOULD BE INCLUDED IN COUNT
    # lon_insill = lon_vals >= sillsea #these are particles in the inner basin or sill any given point in time
    # lon_est = lon_vals >= 0 #these are particles anywhere in the estuary
    # lon_out = (lon_vals <= sillsea) & (lon_vals >= 0) #these are particles in the outer basin at any point in time

    # lon_in_lower = lon_in * lon_in[np.newaxis, 0, :] * z_start_lower #particles in the inner basin that started in the inner basin below sill depth
    # lon_in_upper = lon_in * lon_in[np.newaxis, 0, :] * z_start_upper #particles in the inner basin that started in the inner basin above sill depth


    # lon_start = lon_vals[np.newaxis, 0, :] #starting longitudes of the particles #DON'T ADD NEWAXIS
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
    #whole estuary
    rt_strict_est = np.argmin(lon_est,axis=0).astype('float')
    rt_strict_est = np.where(rt_strict_est==0, tmax+1, rt_strict_est) #replace 0 with tmax+1 (every particle is in the estuary at t=0, these are particles that never leave estuary)
    rt_strict_days_est = rt_strict_est/24 #convert to days
    print('got residence times\n')

    # lon_out_lower = lon_out * lon_out[np.newaxis, 0, :] * z_start_lower #particles in the outer basin that started in the outer basin below sill depth
    # lon_out_upper = lon_out * lon_out[np.newaxis, 0, :] * z_start_upper #particles in the outer basin that started in the outer basin above sill depth

    # print('got lon_out\n')


    # par_in_lower = np.sum(lon_in_lower,axis=1)
    # par_in_upper = np.sum(lon_in_upper,axis=1)
    # par_out_lower = np.sum(lon_out_lower,axis=1)
    # par_out_upper = np.sum(lon_out_upper,axis=1)
    # print('got particle counts\n')
    # lon1 = dsr.lon.where((dsr.lon.sel(Time=0)<sillsea),drop=True) #THESE ARE THE OUTER BASIN PARTICLES
    # print('got lon1\n')
    # sys.stdout.flush()
    #SKIP SILL PARTICLES FOR NOW
    # lon2 = dsr.lon.where((dsr.lon.sel(Time=0)>=sillsea) & (dsr.lon.sel(Time=0)<sillland),drop=True)
    # print('got lon2\n')
    # # sys.stdout.flush()
    # lon3 = dsr.lon.where((dsr.lon.sel(Time=0)>=sillland),drop=True) #THESE ARE THE INNER BASIN PARTICLES
    # print('got lon3\n')
    # sys.stdout.flush()




    #FOR NOW, SKIP PARTICLES INITIALIZED ON SILL, AND ONLY PLOT PARTICLES REMAINING IN ORIGINAL BASIN
    # par1_ocn=(lon1<0).astype(int).sum(dim='Particle')
    # par1_out=((lon1>=0) & (lon1<sillsea)).astype(int).sum(dim='Particle') #THESE ARE THE OUTER BASIN PARTICLES STILL IN OUTER BASIN
    # par1_sill=((lon1>=sillsea) & (lon1<sillland)).astype(int).sum(dim='Particle')
    # par1_in=(lon1>=sillland).astype(int).sum(dim='Particle')
    # print('got par1\n')
    # sys.stdout.flush()

    # par2_ocn=(lon2<0).astype(int).sum(dim='Particle')
    # par2_out=((lon2>=0) & (lon2<sillsea)).astype(int).sum(dim='Particle')
    # par2_sill=((lon2>=sillsea) & (lon2<sillland)).astype(int).sum(dim='Particle')
    # par2_in=(lon2>=sillland).astype(int).sum(dim='Particle')
    # print('got par2\n')
    # sys.stdout.flush()

    # par3_ocn=(lon3<0).astype(int).sum(dim='Particle')
    # par3_out=((lon3>=0) & (lon3<sillsea)).astype(int).sum(dim='Particle')
    # par3_sill=((lon3>=sillsea) & (lon3<sillland)).astype(int).sum(dim='Particle')
    # par3_in=(lon3>=sillland).astype(int).sum(dim='Particle') #THESE ARE THE INNER BASIN PARTICLES STILL IN INNER BASIN
    # print('got par3\n')
    # sys.stdout.flush()

    # lonest = (lon>0)
    # lonest = lonest.astype(int)
    # partest = lonest.sum(dim='Particle')
    # partest_ta = zfun.lowpass(partest.values,f='godin',nanpad=True)[35:-35]
    # tplot = partest.Time.values[35:-35]


    # ax1.plot(par1_ocn.Time/24, zfun.lowpass((par1_ocn/par1_out.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:green', label='Ocean')
    # ax1.plot(par1_out.Time/24, zfun.lowpass((par1_out/par1_out.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:cyan', label='Outer basin')
    # ax1.plot(par1_out.Time/24, zfun.lowpass((par1_out/par1_out.sel(Time=0)).values*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY OUTER PARTICLES REMAINING
    # ax1.plot(par1_sill.Time/24, zfun.lowpass((par1_sill/par1_out.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:purple', label='Sill')
    # ax1.plot(par1_in.Time/24, zfun.lowpass((par1_in/par1_out.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:pink', label='Inner basin')

    #ax.set_ylim(20000,35000)


    # ax2.plot(par2_ocn.Time/24, zfun.lowpass((par2_ocn/par2_sill.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:green', label='Ocean')
    # ax2.plot(par2_out.Time/24, zfun.lowpass((par2_out/par2_sill.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:cyan', label='Outer basin')
    # ax2.plot(par2_sill.Time/24, zfun.lowpass((par2_sill/par2_sill.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:purple', label='Sill')
    # ax2.plot(par2_in.Time/24, zfun.lowpass((par2_in/par2_sill.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:pink', label='Inner basin')


    # ax3.plot(par3_ocn.Time/24, zfun.lowpass((par3_ocn/par3_in.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:green', label='Ocean')
    # ax3.plot(par3_out.Time/24, zfun.lowpass((par3_out/par3_in.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:cyan', label='Outer basin')
    # ax3.plot(par3_sill.Time/24, zfun.lowpass((par3_sill/par3_in.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:purple', label='Sill')
    # ax3.plot(par3_in.Time/24, zfun.lowpass((par3_in/par3_in.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:pink', label='Inner basin')
    # ax3.plot(par3_in.Time/24, zfun.lowpass((par3_in/par3_in.sel(Time=0)).values*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY INNER PARTICLES REMAINING
    

    # axs[0,0].plot(time_hours/24, zfun.lowpass((par_in_lower/par_in_lower[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY INNER PARTICLES REMAINING
    # axs[0,1].plot(time_hours/24, zfun.lowpass((par_in_upper/par_in_upper[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY INNER PARTICLES REMAINING
    # axs[1,0].plot(time_hours/24, zfun.lowpass((par_out_lower/par_out_lower[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY OUTER PARTICLES REMAINING
    # axs[1,1].plot(time_hours/24, zfun.lowpass((par_out_upper/par_out_upper[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY OUTER PARTICLES REMAINING

    # axs[0,0].plot(time_hours/24, zfun.lowpass(par_in_lower, f='godin'), color=linecolor, label=silllenlabel) #TRY WITH TOTAL PARTICLE COUNTS
    # axs[0,1].plot(time_hours/24, zfun.lowpass(par_in_upper, f='godin'), color=linecolor, label=silllenlabel) 
    # axs[1,0].plot(time_hours/24, zfun.lowpass(par_out_lower, f='godin'), color=linecolor, label=silllenlabel) 
    # axs[1,1].plot(time_hours/24, zfun.lowpass(par_out_upper, f='godin'), color=linecolor, label=silllenlabel) 
    
    #bin based on salinity now
    z_bin_edges = np.arange(-200,25,25) #25m increments from -200 to 0, which works with the release layers #try to see if we can plot with positive depths
    # s_bin_edges = np.arange(0,34.5,0.5) #half psu increments from 0 to 34 (max s)
    # s_bin_edges = np.arange(22,34.5,0.5) #half psu increments from 0 to 34 (max s) #could try smaller/bigger bins
    # s_bin_edges = np.arange(22,34.25,0.25) #quarter psu increments from 0 to 34 (max s) #could try smaller/bigger bins
    # s_bin_edges = np.arange(22,35,1) #psu increments from 0 to 34 (max s) #could try smaller/bigger bins
    # lon_bin_edges = lonp[0,:] #this also includes longitudes in the ocean half, but we can crop it out in the plot
    # lon_bin_edges_pos = np.delete(lon_bin_edges,np.where(lon_bin_edges<0))
    # x_bin_centers_km = llxyfun.lon2x((lon_bin_edges[:-1]+lon_bin_edges[1:])/2,0,45)/1000
    # rt_xmean = stats.binned_statistic(lon_start, rt_strict_days, statistic='mean', bins=lon_bin_edges, range=None)

    #in
    rt_zmean_in = stats.binned_statistic(z_start, rt_strict_days_in, statistic=np.nanmean, bins=z_bin_edges, range=None)
    rt_std_in = stats.binned_statistic(z_start, rt_strict_days_in, statistic=np.nanstd, bins=z_bin_edges, range=None) #Try adding standard deviation as error
    z_bin_centers_in = (rt_zmean_in.bin_edges[:-1]+rt_zmean_in.bin_edges[1:])/2 #use the subsampled bin edges for plotting
    # axs[2].fill_between(-z_bin_centers_in,rt_zmean_in.statistic-rt_std_in.statistic,rt_zmean_in.statistic+rt_std_in.statistic,color=linecolor2,alpha=0.3)
    axs[1].plot(-z_bin_centers_in,rt_zmean_in.statistic,color=linecolor,label=silllenlabel)

    #estuary
    rt_zmean_est = stats.binned_statistic(z_start, rt_strict_days_est, statistic=np.nanmean, bins=z_bin_edges, range=None)
    rt_std_est = stats.binned_statistic(z_start, rt_strict_days_est, statistic=np.nanstd, bins=z_bin_edges, range=None) #Try adding standard deviation as error
    z_bin_centers_est = (rt_zmean_est.bin_edges[:-1]+rt_zmean_est.bin_edges[1:])/2 #use the subsampled bin edges for plotting
    # axs[0].fill_between(-z_bin_centers_est,rt_zmean_est.statistic-rt_std_est.statistic,rt_zmean_est.statistic+rt_std_est.statistic,color=linecolor2,alpha=0.3)
    axs[0].plot(-z_bin_centers_est,rt_zmean_est.statistic,color=linecolor,label=silllenlabel)

    # ax.hist(rt_strict_days,bins=[0,10,20,30,40,50,60,70,80,90,100,110,120], density=True, histtype='step',color=linecolor,linewidth=2,label=silllenlabel)

    #could try with total number of particles and/or double axis
    
    # axs[0,0].plot(time_hours/24, (par_out/par_out[0])*100, color=linecolor, label=silllenlabel) #TRY WITH NO FILTERING
    # axs[0,1].plot(time_hours/24, (par_in/par_in[0])*100, color=linecolor, label=silllenlabel)
    print('plotted\n')
    sys.stdout.flush()
    
    #dsg.close()


#plt.show()
#pfun.end_plot()

# #PLOTTING - HISTOGRAMS
# fig, axs = plt.subplots(5,1,sharex=True)
# for j in range(5):
#     hour=j*180
#     axs[j].set_title('t='+str(hour)+'h')
#     axs[j].hist(dsr['lon'].sel(Time=hour),bins=20,alpha=0.5)
#     #axs[j].set_ylim(0, 30)

#plt.show()
# axs[0,0].set_xlabel('Days')
# axs[0,0].set_ylabel('% of particles remaining in inner basin')
# axs[0,0].set_ylabel('Particles remaining in inner basin') #TRY WITH TOTAL PARTICLE COUNT

axs[0].set_ylabel('Average strict residence time [days]')
axs[0].set_title('Whole estuary')
axs[1].set_title('Inner basin')

axs[0].grid(True)
axs[1].grid(True)
# axs[0,0].set_ylim(0,100)
# axs[0,0].set_ylim(0,par_in_lower[0])
axs[0].set_xlabel('Release depth [m]')
axs[1].set_xlabel('Release depth [m]')
axs[0].legend(loc='upper left')
axs[1].legend(loc='upper left')
axs[0].set_box_aspect(1)
axs[1].set_box_aspect(1)
axs[0].text(.02, .98, 'A', horizontalalignment='left', verticalalignment='top', transform=axs[0].transAxes, fontsize=14, fontweight='bold')
axs[1].text(.02, .98, 'B', horizontalalignment='left', verticalalignment='top', transform=axs[1].transAxes, fontsize=14, fontweight='bold')

axs[0].set_xlim(0,200) #x axis is the depth
axs[0].set_ylim(0,80)
# axs[0].invert_xaxis()

axs[1].set_xlim(0,200)
axs[1].set_ylim(0,80)
# axs[2].invert_xaxis()


# axs[0,1].set_xlabel('Days')
# # axs[0,1].set_ylabel('% of particles')
# axs[0,1].set_title('Released in inner basin above sill height')
# axs[0,1].grid(True)
# axs[0,1].set_xlim(0,120)
# # axs[0,1].set_ylim(0,100)
# axs[0,1].set_ylim(0,par_in_upper[0])


# axs[1,0].set_ylabel('% of particles remaining in outer basin')
# axs[1,0].set_ylabel('Particles remaining in outer basin') #TRY WITH TOTAL PARTICLE COUNT
# axs[1,0].set_title('Released in outer basin below sill height')
# axs[1,0].grid(True)
# axs[1,0].set_xlim(0,120)
# # axs[1,0].set_ylim(0,100)
# axs[1,0].set_ylim(0,par_out_lower[0])

# axs[1,1].set_xlabel('Days')
# # axs[1,1].set_ylabel('% of particles')
# axs[1,1].set_title('Released in outer basin above sill height')
# axs[1,1].grid(True)
# axs[1,1].set_xlim(0,120)
# # axs[1,1].set_ylim(0,100)
# axs[1,1].set_ylim(0,par_out_upper[0])


# ax2.set_xlabel('Days')
# #ax1.set_ylabel('Number of particles')
# ax2.set_title('Particles released on sill')
# #ax2.legend(loc='best')
# ax2.grid(True)
# ax2.set_xlim(0,120)
# ax2.set_ylim(0,100)

# ax3.set_xlabel('Days')
# #ax3.set_ylabel('Number of particles')
# ax3.set_title('Particles released in inner basin (unfiltered)')
# #ax3.legend(loc='best')
# ax3.grid(True)
# ax3.set_xlim(0,120)
# ax3.set_ylim(0,100)
# ax3.legend(loc='upper right')


fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_tracker2_zdist_regions2.png'
plt.savefig(fn_fig)
plt.close()
#plt.show()

