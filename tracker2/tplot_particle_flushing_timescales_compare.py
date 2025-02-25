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
import tef_fun
import datetime
import pandas as pd

plt.close('all')
# fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(20,6))
# fig, axs = plt.subplots(2,2,figsize=(15,8), gridspec_kw={'height_ratios': [6,1]})

efold_est_days=np.zeros(5)
efold_inner_days=np.zeros(5)
efold_noreentry_est_days=np.zeros(5)
efold_noreentry_inner_days=np.zeros(5)
residence_est_days=np.zeros(5)
residence_inner_days=np.zeros(5)
exposure_est_days=np.zeros(5)
exposure_inner_days=np.zeros(5)

silllens_plot=[5,10,20,40,80]

for i in range(5):
    # Choose an experiment and release to plot.
    # in_dir0 = Ldir['LOo'] / 'tracks'
    # exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    #     itext='** Choose experiment from list **', last=False)
    # rel = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    #     itext='** Choose item from list **', last=False)

    #loop through five models
    if i==0:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(45e3,0,45)
        linecolor = 'tab:red'
        linecolor2 = plt.cm.tab20(7)
        silllenlabel = '5km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill5km_t0_xa0/sill5kmest_3d/release_2020.09.01.nc'
        tef_5km_fn = Ldir['LOo'] / 'extract/sill5km_t0_xa0/tef2/bulk_hourly_2020.09.01_2020.12.31' #this is for the Qprism timekeeper and shading
    elif i==1:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(50e3,0,45)
        linecolor = 'tab:orange'
        linecolor2 = plt.cm.tab20(3)
        silllenlabel = '10km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill10km_t2_xa0/sill10kmest_3d/release_2020.09.01.nc'
    elif i==2:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(60e3,0,45)
        linecolor = 'tab:green'
        linecolor2 = plt.cm.tab20(5)
        silllenlabel = '20km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill20kmdeep_t2_xa0/sill20kmdeepest_3d/release_2020.09.01.nc'
    elif i==3:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(80e3,0,45)
        linecolor = 'tab:blue'
        linecolor2 = plt.cm.tab20(1)
        silllenlabel = '40km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill40km_t2_xa0/sill40kmest_3d/release_2020.09.01.nc'
    elif i==4:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(120e3,0,45)
        linecolor = 'tab:purple'
        linecolor2 = plt.cm.tab20(9)
        silllenlabel = '80km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill80km_t2_xa0/sill80kmest_3d/release_2020.09.01.nc'
    
    print(silllenlabel+'\n')

    #dataset
    dsr = xr.open_dataset(fn, decode_times=False)
    NT, NP = dsr.lon.shape
    #longitude data
    lon_vals = dsr.lon.values #longitudes of the particles
    #time data
    time_hours = dsr.Time.values
    dsr.close()
    print('got lon_vals and time\n')

    #longitude data
    lon_start = lon_vals[np.newaxis, 0, :] #starting longitudes of the particles
    lon_start_in = lon_start >= sillland #boolean array for particles starting in the inner basin
    # lon_start_insill = lon_start >= sillsea #boolean array for particles starting in the inner basin and sill

    lon_in = lon_vals >= sillland #boolean array of particles in the inner basin over time
    # lon_insill = lon_vals >= sillsea #boolean array of particles in the inner basin over time
    lon_est = lon_vals >= 0 #boolean array of particles in the whole estuary over time

    #this includes returning particles
    start_in_stay_in = lon_in * lon_start_in #particles in the inner basin that started in the inner basin
    # start_insill_stay_insill = lon_insill * lon_start_insill #particles in the inner basin and sill that started in the inner basin or sill #could change this to use lon_start_in?
    start_est_stay_est = lon_est

    #to count particles that have never left, find the time of first exit
    tmax = time_hours[-1]
    #inner basin
    rt_strict_in = np.argmin(lon_in,axis=0).astype('float') #first time the particle is outside the inner basin
    rt_strict_in = np.where(rt_strict_in==0, tmax+1, rt_strict_in) #replace 0 with tmax+1, this sets the particles that never leave or ones that were released outside to tmax+1
    rt_strict_in = rt_strict_in * lon_start_in #this resets the particles that are not released in the inner basin to zero (necessary?)
    rt_ind_in = rt_strict_in.astype(int)
    exit_mask_in = np.zeros(lon_vals.shape, bool)
    exit_mask_in[rt_ind_in[rt_ind_in<tmax],np.flatnonzero(rt_ind_in<tmax)] = True #this is a mask with the first time a particle exits set as True
    strict_mask_in = np.logical_not(np.cumsum(exit_mask_in, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False
    #inner basin + sill
    # rt_strict_insill = np.argmin(lon_insill,axis=0).astype('float')
    # rt_strict_insill = np.where(rt_strict_insill==0, tmax+1, rt_strict_insill) #replace 0 with tmax+1, this sets the particles that never leave or ones that were released outside to tmax+1
    # rt_strict_insill = rt_strict_insill * lon_start_insill #this resets the particles that are not released in the inner basin or sill to zero (necessary?)
    # rt_ind_insill = rt_strict_insill.astype(int)
    # exit_mask_insill = np.zeros(lon_vals.shape, bool)
    # exit_mask_insill[rt_ind_insill[rt_ind_insill<tmax],np.flatnonzero(rt_ind_insill<tmax)] = True #this is a mask with the first time a particle exits set as True
    # strict_mask_insill = np.logical_not(np.cumsum(exit_mask_insill, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False
    #whole estuary
    rt_strict_est = np.argmin(lon_est,axis=0).astype('float')
    rt_strict_est = np.where(rt_strict_est==0, tmax+1, rt_strict_est) #replace 0 with tmax+1 (every particle is in the estuary at t=0, these are particles that never leave estuary)
    rt_ind_est = rt_strict_est.astype(int)
    exit_mask_est = np.zeros(lon_vals.shape, bool)
    exit_mask_est[rt_ind_est[rt_ind_est<tmax],np.flatnonzero(rt_ind_est<tmax)] = True #this is a mask with the first time a particle exits set as True
    strict_mask_est = np.logical_not(np.cumsum(exit_mask_est, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False

    start_in_stay_in_strict = start_in_stay_in*strict_mask_in
    # start_insill_stay_insill_strict = start_insill_stay_insill*strict_mask_insill
    start_est_stay_est_strict = start_est_stay_est*strict_mask_est

    print('got data and boolean arrays\n')

    #get particle counts
    count_in = np.sum(start_in_stay_in,axis=1) #total particles in the inner basin that started in the inner basin
    # count_insill = np.sum(start_insill_stay_insill,axis=1) #total particles in the inner basin and sill that started in the inner basin or sill
    count_est = np.sum(start_est_stay_est,axis=1) #total particles in the estuary

    count_in_strict = np.sum(start_in_stay_in_strict,axis=1) #particles that have never left the inner basin
    # count_insill_strict = np.sum(start_insill_stay_insill_strict,axis=1) #particles that have never left the inner basin + sill
    count_est_strict = np.sum(start_est_stay_est_strict,axis=1) #particles that have never left the estuary
    print('got particle counts\n')
    sys.stdout.flush()

    #tidally average and scale the data
    pad=36 #could use 35 but 36 is nice bc it gives exactly 1.5 days removed on each end

    count_in_ta = zfun.lowpass(count_in, f='godin')[pad:-pad+1] 
    count_est_ta = zfun.lowpass(count_est, f='godin')[pad:-pad+1]
    count_in_strict_ta = zfun.lowpass(count_in_strict, f='godin')[pad:-pad+1] 
    count_est_strict_ta = zfun.lowpass(count_est_strict, f='godin')[pad:-pad+1]

    frac_in_ta = zfun.lowpass((count_in/count_in[0]), f='godin')[pad:-pad+1] 
    frac_est_ta = zfun.lowpass((count_est/count_est[0]), f='godin')[pad:-pad+1]
    frac_in_strict_ta = zfun.lowpass((count_in_strict/count_in_strict[0]), f='godin')[pad:-pad+1] 
    frac_est_strict_ta = zfun.lowpass((count_est_strict/count_est_strict[0]), f='godin')[pad:-pad+1]

    par_scale_in = count_in_ta[0]
    par_scale_est = count_est_ta[0]
    par_scale_in_strict = count_in_strict_ta[0]
    par_scale_est_strict = count_est_strict_ta[0]

    t_ta = time_hours[pad:-pad+1]
    t_scale = t_ta[-1]
    t_frac = t_ta / t_scale

    #exponential fit with two parameters
    def func2(x, a, b):
        return a * np.exp(-x/b) #change b to be time units

    p0_guess=(1,1) #because the data has been scaled, guess 1,1 as the parameters a and b
    popt_in, pcov_in = curve_fit(func2, t_frac, frac_in_ta, p0=p0_guess)
    popt_est, pcov_est = curve_fit(func2, t_frac, frac_est_ta, p0=p0_guess)
    popt_in_strict, pcov_in_strict = curve_fit(func2, t_frac, frac_in_strict_ta, p0=p0_guess)
    popt_est_strict, pcov_est_strict = curve_fit(func2, t_frac, frac_est_strict_ta, p0=p0_guess)

    T_e_in_fit_days = (popt_in[1])*t_scale/24 #T_e from two parameter fit in days
    T_e_est_fit_days = (popt_est[1])*t_scale/24 #T_e from two parameter fit in days
    T_e_in_strict_fit_days = (popt_in_strict[1])*t_scale/24 #T_e from two parameter fit in days
    T_e_est_strict_fit_days = (popt_est_strict[1])*t_scale/24 #T_e from two parameter fit in days      
    A_in_fit = popt_in[0]
    A_est_fit = popt_est[0]
    A_in_strict_fit = popt_in_strict[0]
    A_est_strict_fit = popt_est_strict[0]

    print('\nT_e_in from two-param fit:\n')
    print(T_e_in_fit_days)
    print('\nT_e_est from two-param fit:\n')
    print(T_e_est_fit_days)
    print('\nT_e_in_strict from two-param fit:\n')
    print(T_e_in_strict_fit_days)
    print('\nT_e_est_strict from two-param fit:\n')
    print(T_e_est_strict_fit_days)

    efold_est_days[i]=T_e_est_fit_days
    efold_inner_days[i]=T_e_in_fit_days
    efold_noreentry_est_days[i]=T_e_est_strict_fit_days
    efold_noreentry_inner_days[i]=T_e_in_strict_fit_days

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
    print('got exposure times\n')
    print('average estuary exposure time [days]')
    print(np.nanmean(exposuret_days_est))
    print('average inner basin exposure time [days]')
    print(np.nanmean(exposuret_days_in))
    residence_est_days[i]=np.nanmean(rt_strict_days_est)
    residence_inner_days[i]=np.nanmean(rt_strict_days_in)
    exposure_est_days[i]=np.nanmean(exposuret_days_est)
    exposure_inner_days[i]=np.nanmean(exposuret_days_in)

#Flushing times
sect_est ='a1'
sect_in='b5'
flushing_est=np.zeros(5)
flushing_inner=np.zeros(5)
flushing_est_days=np.zeros(5)
flushing_inner_days=np.zeros(5)
silllens=['5km', '10km', '20km', '40km', '80km']
gridnames = ['sill5km', 'sill10km', 'sill20kmdeep', 'sill40km', 'sill80km']
gctags=['sill5km_c0', 'sill10km_c0', 'sill20kmdeep_c0', 'sill40km_c0', 'sill80km_c0']
gtagexs=['sill5km_t0_xa0', 'sill10km_t2_xa0', 'sill20kmdeep_t2_xa0', 'sill40km_t2_xa0', 'sill80km_t2_xa0']
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'
Ldir['ds0']='2020.09.01'
Ldir['ds1']='2020.12.31'
# take a subset of the data so that we are averaging over an integer number of spring neap cycles at the end
# use 7 spring neap cycles, starting at index 257 and going to 2741 - these are the peaks in the 5km Qprism but similar for the other models
start_avg_ind = 257
end_avg_ind = 2741
start_avg_dt = '2020-09-13 05:00:00'
end_avg_dt = '2020-12-25 16:00:00' #this is tef_df.index[2740], the time indexing doesn't seem to be open half interval
for i in range(5):
# for i in range(len(gctags)-1):
# for i in [0,2]:
    #model and extraction info
    print(silllens[i])
    gctag=gctags[i]
    gtagex=gtagexs[i]
    gridname=gridnames[i]
    in_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('bulk_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
    seg_est_inner_fn = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('segments_' + Ldir['ds0'] + '_' + Ldir['ds1'] + '_' + gridname + '_cei_rivA1.nc')
    
    #Whole estuary
    sect_name = sect_est #a1
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)         
    # adjust units
    # tef_df['Q_p'] = tef_df['q_p']/1000
    # tef_df['Q_m'] = tef_df['q_m']/1000
    # tef_df['Q_prism']=tef_df['qprism']/1000
    # get variables for efflux-reflux calculation
    Qout_est = -tef_df['q_m'][start_avg_ind:end_avg_ind].mean() #keep original units for comparing with V
    sout_est = tef_df['salt_m'][start_avg_ind:end_avg_ind].mean()

    #Section b5
    sect_name = sect_in #b5
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
    # adjust units
    # tef_df['Q_p'] = tef_df['q_p']/1000
    # tef_df['Q_m'] = tef_df['q_m']/1000
    # tef_df['Q_prism']=tef_df['qprism']/1000
    # get variables for efflux-reflux calculation
    Qout_inner = -tef_df['q_m'][start_avg_ind:end_avg_ind].mean()
    sout_inner = tef_df['salt_m'][start_avg_ind:end_avg_ind].mean()

    #Get the volume of the estuary and inner basin
    #Might be able to get this from the segments?
    grid_dir = Ldir['data'] / 'grids' / gridname
    #get grid info
    g = xr.open_dataset(grid_dir / 'grid.nc')
    h = g.h.values
    mask_rho=g.mask_rho.values
    h[g.mask_rho.values==0] = np.nan
    lon_rho = g.lon_rho.values
    lat_rho = g.lat_rho.values
    lon_u = g.lon_u.values
    DX = 1/g.pm.values
    DY = 1/g.pn.values
    DA = DX*DY
    #find the volume between the two longitudes
    #get the actual longitude of the sections
    sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
    sect_df = pd.read_pickle(sect_df_fn)
    #get the a1 part of the sect_df
    sect_df_a1=sect_df[sect_df.sn=='a1']
    #get the index of the u values
    a1_ind_u=sect_df_a1.i.values[0]
    #get the longitude of the section from lon_u
    a1_lon_u=lon_u[0,a1_ind_u]
    #do the same for b5
    sect_df_b5=sect_df[sect_df.sn=='b5']
    b5_ind_u=sect_df_b5.i.values[0]
    b5_lon_u=lon_u[0,b5_ind_u]
    #to get estuary volume mask the h array ouside of a1
    h_est = g.h.values
    h_est[mask_rho==0] = np.nan #mask the land
    h_est[lon_rho<a1_lon_u] = np.nan #mask seaward of inner basin
    V_est = np.nansum(DA*h_est)
    #to get inner basin volume mask the h array ouside of b5
    h_inner = g.h.values
    h_inner[mask_rho==0] = np.nan #mask the land
    h_inner[lon_rho<b5_lon_u] = np.nan #mask seaward of inner basin
    V_inner = np.nansum(DA*h_inner)

    #Calculate the flushing times
    flushing_est[i] = V_est/Qout_est
    flushing_inner[i] = V_inner/Qout_inner
    flushing_est_days[i] = flushing_est[i]/(24*3600)
    flushing_inner_days[i] = flushing_inner[i]/(24*3600)
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

# axs[0,0].set_ylabel('% of particles')
# axs[0,0].set_title('Particles released in outer basin')
# # ax1.set_title('Particles released in outer basin (unfiltered)')
# #ax1.legend(loc='best')
# axs[0,0].grid(True)
# axs[0,0].set_xlim(0,120)
# axs[0,0].set_ylim(0,100)

# #legend
# handles, labels = axs[0,0].get_legend_handles_labels()
# handles_reorder = np.concatenate((handles[::2],handles[1::2]),axis=0)
# labels_reorder = np.concatenate((labels[::2],labels[1::2]),axis=0)
# axs[0,0].legend(handles_reorder,labels_reorder,loc='upper right',ncol=2)

# # ax2.set_xlabel('Days')
# # #ax1.set_ylabel('Number of particles')
# # ax2.set_title('Particles released on sill')
# # #ax2.legend(loc='best')
# # ax2.grid(True)
# # ax2.set_xlim(0,120)
# # ax2.set_ylim(0,100)

# # axs[0,1].set_xlabel('Days')
# #ax3.set_ylabel('Number of particles')
# axs[0,1].set_title('Particles released in inner basin')
# # ax3.set_title('Particles released in inner basin (unfiltered)')
# #ax3.legend(loc='best')
# axs[0,1].grid(True)
# axs[0,1].set_xlim(0,120)
# axs[0,1].set_ylim(0,100)

# axs[1,0].set_xlim(0,120)
# axs[1,1].set_xlim(0,120)
# axs[1,0].set_ylim(20,80)
# axs[1,1].set_ylim(20,80)
# axs[1,0].set_xlabel('Days')
# axs[1,1].set_xlabel('Days')

# fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_tracker2_fitandcalc_sn.png'
# plt.savefig(fn_fig)
# plt.close()
# #plt.show()

#update plot style for thesis
# pfun.start_plot(fs=14)
fig, [ax1,ax2] = plt.subplots(1,2,figsize=(12,6))
# fig, [ax1,ax2] = plt.subplots(1,2,figsize=(8,5))
ax1.plot(silllens_plot,exposure_est_days,c='tab:brown',marker='o',label='Average exposure time')
ax1.plot(silllens_plot,residence_est_days,c=plt.cm.tab20b(14),marker='o',label='Average strict residence time')
ax1.plot(silllens_plot,efold_est_days,c=plt.cm.tab20b(10),marker='o',label='e-folding time')
ax1.plot(silllens_plot,efold_noreentry_est_days,c=plt.cm.tab20b(4),marker='o',label='e-folding time (no reentry)')
ax1.plot(silllens_plot,flushing_est_days,c='teal',marker='o',label='Flushing time')
ax2.plot(silllens_plot,exposure_inner_days,c='tab:brown',marker='o',label='Average exposure time')
ax2.plot(silllens_plot,residence_inner_days,c=plt.cm.tab20b(14),marker='o',label='Average strict residence time')
ax2.plot(silllens_plot,efold_inner_days,c=plt.cm.tab20b(10),marker='o',label='e-folding time')
ax2.plot(silllens_plot,efold_noreentry_inner_days,c=plt.cm.tab20b(4),marker='o',label='e-folding time (no reentry)')
ax2.plot(silllens_plot,flushing_inner_days,c='teal',marker='o',label='Flushing time')
ax1.set_xlabel('Sill length [km]')
ax2.set_xlabel('Sill length [km]')
ax1.set_ylabel('Time scale [days]')
ax1.set_title('Whole estuary')
ax2.set_title('Inner basin')
ax1.grid(True)
ax2.grid(True)
ax1.legend(loc='lower right')
# ax2.legend(loc='lower right')
ax1.set_xlim(0,85)
ax2.set_xlim(0,85)
ax1.set_ylim(0,180) #set these once we see the results
ax2.set_ylim(0,180)
ax1.set_box_aspect(1)
ax2.set_box_aspect(1)
ax1.text(.02, .98, 'A', horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes, fontsize=14, fontweight='bold')
ax2.text(.02, .98, 'B', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes, fontsize=14, fontweight='bold')
fn_fig = Ldir['LOo'] / 'plots' / 'particle_flushing_timescales_compare.png'
plt.savefig(fn_fig)
plt.close()
