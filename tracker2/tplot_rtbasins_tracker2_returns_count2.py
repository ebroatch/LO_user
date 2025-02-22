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

plt.close('all')
# fig, axs = plt.subplots(1,5,figsize=(25,5),sharey=True)#,gridspec_kw={'height_ratios': [6,1]})
fig, axs = plt.subplots(1,5,figsize=(18,6),sharey=True,layout='constrained')#,gridspec_kw={'height_ratios': [6,1]})
# fig, axs = plt.subplots(2,2,figsize=(15,15))

for i in range(5):
    # Choose an experiment and release to plot.
    # in_dir0 = Ldir['LOo'] / 'tracks'
    # exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    #     itext='** Choose experiment from list **', last=False)
    # rel = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    #     itext='** Choose item from list **', last=False)


    #for now just use 5,20,80 - can add 10 and 40 if it runs fast
    if i==0:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(45e3,0,45)
        linecolor = 'tab:red'
        silllenlabel = '5km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill5km_t0_xa0/sill5kmest_3d/release_2020.09.01.nc'
        tef_5km_fn = Ldir['LOo'] / 'extract/sill5km_t0_xa0/tef2/bulk_hourly_2020.09.01_2020.12.31' #this is for the Qprism timekeeper and shading
    elif i==1:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(50e3,0,45)
        linecolor = 'tab:orange'
        silllenlabel = '10km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill10km_t2_xa0/sill10kmest_3d/release_2020.09.01.nc'
    elif i==2:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(60e3,0,45)
        linecolor = 'tab:green'
        silllenlabel = '20km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill20kmdeep_t2_xa0/sill20kmdeepest_3d/release_2020.09.01.nc'
    elif i==3:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(80e3,0,45)
        linecolor = 'tab:blue'
        silllenlabel = '40km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill40km_t2_xa0/sill40kmest_3d/release_2020.09.01.nc'
    elif i==4:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(120e3,0,45)
        linecolor = 'tab:purple'
        silllenlabel = '80km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill80km_t2_xa0/sill80kmest_3d/release_2020.09.01.nc'
    
    print(silllenlabel+'\n')

    # get Datasets
    #fn = in_dir0 / exp_name / rel
    #fng = in_dir0 / exp_name / 'grid.nc'
    dsr = xr.open_dataset(fn, decode_times=False)
    #dsg = xr.open_dataset(fng)

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

    #longitude data
    lon_vals = dsr.lon.values #longitudes of the particles
    # z_vals = dsr.z.values #depths of the particles
    
    time_hours = dsr.Time.values
    dsr.close()
    print('got lon_vals and time\n')

    lon_start = lon_vals[np.newaxis, 0, :] #starting longitudes of the particles
    lon_start_in = lon_start >= sillland #boolean array for particles starting in the inner basin
    lon_start_insill = lon_start >= sillsea #boolean array for particles starting in the inner basin and sill
    # lon_start_out = (lon_start <= sillsea) & (lon_start >= 0) #boolean array for particles starting in the outer basin

    # z_start = z_vals[np.newaxis, 0, :] #starting depth of the particles
    # z_start_low = z_start < -50 #boolean array for particles starting below sill depth
    # z_start_up = z_start >= -50 #boolean array for particles starting above sill depth

    lon_in = lon_vals >= sillland #boolean array of particles in the inner basin over time
    lon_insill = lon_vals >= sillsea #boolean array of particles in the inner basin over time
    # lon_est = lon_vals >= 0 #boolean array of particles in the whole estuary over time
    # lon_out = (lon_vals <= sillsea) & (lon_vals >= 0) #boolean array of particles in the outer basin over time

    # z_low = z_vals < -50 #boolean array of particles below sill depth over time
    # z_up = z_vals >= -50 #boolean array of particles above sill depth over time

    #this includes returning particles
    start_in_stay_in = lon_in * lon_start_in #particles in the inner basin that started in the inner basin
    start_insill_stay_insill = lon_insill * lon_start_insill #particles in the inner basin and sill that started in the inner basin or sill #could change this to use lon_start_in?
    # start_est_stay_est = lon_est
    # start_out_stay_out = lon_out * lon_start_out #particles in the outer basin that started in the outer basin

    #count number of returns to inner basin
    ret_exit = np.zeros(lon_in.shape) 
    ret_exit[1:,:] = np.diff(lon_in.astype(int), axis=0) #-1 for exit, +1 for return !!!important to cast lon_in as int or this won't work, it will count exits and returns as 1!!!
    ret = np.where(ret_exit==1, 1, 0) #only keep the returns as +1 for the first hour that the particle is back inside the domain
    ret_total = np.cumsum(ret, axis=0) #the number of returns a particle has made
    ret_in_total = np.where(lon_in==0,np.nan,ret_total) #remove the particles that are currently outside
    ret_in_total = np.where(lon_start_in==0,np.nan,ret_in_total) #remove the particles that started outside

    ret_0 = ret_in_total==0
    ret_1 = ret_in_total==1
    ret_2 = ret_in_total==2
    ret_3 = ret_in_total==3
    ret_4plus = ret_in_total>=4

    print('got data and boolean arrays\n')
    
    ret_count_0 = np.nansum(ret_0,axis=1)
    ret_count_1 = np.nansum(ret_1,axis=1)
    ret_count_2 = np.nansum(ret_2,axis=1)
    ret_count_3 = np.nansum(ret_3,axis=1)
    ret_count_4plus = np.nansum(ret_4plus,axis=1)

    ret_count_0_ta = zfun.lowpass(ret_count_0, f='godin')
    ret_count_1_ta = zfun.lowpass(ret_count_1, f='godin')
    ret_count_2_ta = zfun.lowpass(ret_count_2, f='godin')
    ret_count_3_ta = zfun.lowpass(ret_count_3, f='godin')
    ret_count_4plus_ta = zfun.lowpass(ret_count_4plus, f='godin')



    # #to count particles that have never left, find the time of first exit
    # tmax = time_hours[-1]
    # #inner basin
    # rt_strict_in = np.argmin(lon_in,axis=0).astype('float') #first time the particle is outside the inner basin
    # rt_strict_in = np.where(rt_strict_in==0, tmax+1, rt_strict_in) #replace 0 with tmax+1, this sets the particles that never leave or ones that were released outside to tmax+1
    # rt_strict_in = rt_strict_in * lon_start_in #this resets the particles that are not released in the inner basin to zero (necessary?)
    # rt_ind_in = rt_strict_in.astype(int)
    # exit_mask_in = np.zeros(lon_vals.shape, bool)
    # exit_mask_in[rt_ind_in[rt_ind_in<tmax],np.flatnonzero(rt_ind_in<tmax)] = True #this is a mask with the first time a particle exits set as True
    # strict_mask_in = np.logical_not(np.cumsum(exit_mask_in, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False
    # #inner basin + sill
    # rt_strict_insill = np.argmin(lon_insill,axis=0).astype('float')
    # rt_strict_insill = np.where(rt_strict_insill==0, tmax+1, rt_strict_insill) #replace 0 with tmax+1, this sets the particles that never leave or ones that were released outside to tmax+1
    # rt_strict_insill = rt_strict_insill * lon_start_insill #this resets the particles that are not released in the inner basin or sill to zero (necessary?)
    # rt_ind_insill = rt_strict_insill.astype(int)
    # exit_mask_insill = np.zeros(lon_vals.shape, bool)
    # exit_mask_insill[rt_ind_insill[rt_ind_insill<tmax],np.flatnonzero(rt_ind_insill<tmax)] = True #this is a mask with the first time a particle exits set as True
    # strict_mask_insill = np.logical_not(np.cumsum(exit_mask_insill, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False
    # #whole estuary
    # rt_strict_est = np.argmin(lon_est,axis=0).astype('float')
    # rt_strict_est = np.where(rt_strict_est==0, tmax+1, rt_strict_est) #replace 0 with tmax+1 (every particle is in the estuary at t=0, these are particles that never leave estuary)
    # rt_ind_est = rt_strict_est.astype(int)
    # exit_mask_est = np.zeros(lon_vals.shape, bool)
    # exit_mask_est[rt_ind_est[rt_ind_est<tmax],np.flatnonzero(rt_ind_est<tmax)] = True #this is a mask with the first time a particle exits set as True
    # strict_mask_est = np.logical_not(np.cumsum(exit_mask_est, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False


    # rt_strict = np.argmax(lon_vals<0,axis=0).astype('float')
    # rt_strict = np.where(rt_strict==0, tmax+1, rt_strict) #replace 0 with tmax+1 (every particle is in the estuary at t=0, these are particles that never leave estuary)
    # rt_ind = rt_strict.astype(int)
    # exit_mask = np.zeros(lon_vals.shape, bool)
    # exit_mask[np.flatnonzero(rt_ind<tmax), rt_ind[rt_ind<tmax]] = True #this is a mask with the first time a particle exits set as True
    # strict_mask = np.logical_not(np.cumsum(exit_mask, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False

    # start_in_stay_in_strict = start_in_stay_in*strict_mask_in
    # start_insill_stay_insill_strict = start_insill_stay_insill*strict_mask_insill
    # start_est_stay_est_strict = start_est_stay_est*strict_mask_est
    # start_out_stay_out_strict = start_in_stay_in*strict_mask #OH NO NEED TO CHANGE THE LOGIC HERE #THIS WILL ONLY WORK FOR THE WHOLE ESTUARY SINCE THAT IS HOW THE RT IS DEFINED!!!

    # start_inlow_stay_inlow = lon_in * z_low * lon_start_in * z_start_low #particles that started in the inner basin below sill depth and stayed there
    # start_inup_stay_inup = lon_in * z_up * lon_start_in * z_start_up #particles that started in the inner basin above sill depth and stayed there
    # start_outlow_stay_outlow = lon_out * z_low * lon_start_out * z_start_low #particles that started in the outer basin below sill depth and stayed there
    # start_outup_stay_outup = lon_out * z_up * lon_start_out * z_start_up #particles that started in the outer basin above sill depth and stayed there

    # lon_in_lower = lon_in * lon_in[np.newaxis, 0, :] * z_start_lower #particles in the inner basin that started in the inner basin below sill depth
    # lon_in_upper = lon_in * lon_in[np.newaxis, 0, :] * z_start_upper #particles in the inner basin that started in the inner basin above sill depth

    # lon_out_lower = lon_out * lon_out[np.newaxis, 0, :] * z_start_lower #particles in the outer basin that started in the outer basin below sill depth
    # lon_out_upper = lon_out * lon_out[np.newaxis, 0, :] * z_start_upper #particles in the outer basin that started in the outer basin above sill depth


    # par_in_lower = np.sum(lon_in_lower,axis=1)
    # par_in_upper = np.sum(lon_in_upper,axis=1)
    # par_out_lower = np.sum(lon_out_lower,axis=1)
    # par_out_upper = np.sum(lon_out_upper,axis=1)

    # count_start_inlow_stay_in = np.sum(start_inlow_stay_in,axis=1) #total particles in the inner basin that started in the inner basin below sill depth
    # count_start_inup_stay_in = np.sum(start_inup_stay_in,axis=1) #total particles in the inner basin that started in the inner basin above sill depth
    # count_start_outlow_stay_out = np.sum(start_outlow_stay_out,axis=1) #total particles in the outer basin that started in the outer basin below sill depth
    # count_start_outup_stay_out = np.sum(start_outup_stay_out,axis=1) #total particles in the outer basin that started in the outer basin above sill depth

    # count_start_inlow_stay_inlow = np.sum(start_inlow_stay_inlow,axis=1) #total particles that started in the inner basin below sill depth and stayed there
    # count_start_inup_stay_inup = np.sum(start_inup_stay_inup,axis=1) #total particles that started in the inner basin above sill depth and stayed there
    # count_start_outlow_stay_outlow = np.sum(start_outlow_stay_outlow,axis=1) #total particles that started in the outer basin below sill depth and stayed there
    # count_start_outup_stay_outup = np.sum(start_outup_stay_outup,axis=1) #total particles that started in the outer basin above sill depth and stayed there

    # count_in = np.sum(start_in_stay_in,axis=1) #total particles in the inner basin that started in the inner basin
    # count_insill = np.sum(start_insill_stay_insill,axis=1) #total particles in the inner basin and sill that started in the inner basin or sill
    # count_est = np.sum(start_est_stay_est,axis=1) #total particles in the estuary

    # count_in_strict = np.sum(start_in_stay_in_strict,axis=1) #particles that have never left the inner basin
    # count_insill_strict = np.sum(start_insill_stay_insill_strict,axis=1) #particles that have never left the inner basin + sill
    # count_est_strict = np.sum(start_est_stay_est_strict,axis=1) #particles that have never left the estuary


    print('got particle counts\n')
    # lon1 = dsr.lon.where((dsr.lon.sel(Time=0)<sillsea),drop=True) #THESE ARE THE OUTER BASIN PARTICLES
    # print('got lon1\n')
    sys.stdout.flush()
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

    # axs[0,0].plot(time_hours/24, zfun.lowpass(count_start_inlow_stay_in, f='godin'), color=linecolor, label=silllenlabel) #TRY WITH TOTAL PARTICLE COUNTS AND ADD STRICT LAYER SORTING
    # axs[0,0].plot(time_hours/24, zfun.lowpass(count_start_inlow_stay_inlow, f='godin'), '--', color=linecolor, label=silllenlabel)
    # axs[0,1].plot(time_hours/24, zfun.lowpass(count_start_inup_stay_in, f='godin'), color=linecolor, label=silllenlabel) 
    # axs[0,1].plot(time_hours/24, zfun.lowpass(count_start_inup_stay_inup, f='godin'), '--', color=linecolor, label=silllenlabel) 
    # axs[1,0].plot(time_hours/24, zfun.lowpass(count_start_outlow_stay_out, f='godin'), color=linecolor, label=silllenlabel) 
    # axs[1,0].plot(time_hours/24, zfun.lowpass(count_start_outlow_stay_outlow, f='godin'), '--', color=linecolor, label=silllenlabel) 
    # axs[1,1].plot(time_hours/24, zfun.lowpass(count_start_outup_stay_out, f='godin'), color=linecolor, label=silllenlabel+' in basin') 
    # axs[1,1].plot(time_hours/24, zfun.lowpass(count_start_outup_stay_outup, f='godin'), '--', color=linecolor, label=silllenlabel+' in layer')

    # #add shading #SKIP SHADING FOR NOW, COULD ADD BACK
    # if i==0:
    #     sect_name='b3'
    #     pad=36
    #     tef_df, vn_list, vec_list = tef_fun.get_two_layer(tef_5km_fn, sect_name)
    #     tef_df['Q_prism']=tef_df['qprism']/1000
    #     Qprism = tef_df['Q_prism'].loc['2020-09-04':'2020-12-28'] #cut off extra pad because qprism uses two godin filters
    #     ot=tef_df.loc['2020-09-04':'2020-12-28'].index
    #     ot_hours_delta = (((ot - datetime.datetime(2020,9,1,0,0,0)).total_seconds())/3600).to_numpy()
    #     axs[0,0].set_ylim(0,42000)
    #     axs[0,1].set_ylim(0,21000)
    #     axs[0,2].set_ylim(0,21000)
    #     axs[1,0].plot(ot_hours_delta/24,Qprism.to_numpy(), color='tab:gray', linewidth=2) #cut off the weird ends
    #     axs[1,1].plot(ot_hours_delta/24,Qprism.to_numpy(), color='tab:gray', linewidth=2)
    #     axs[1,2].plot(ot_hours_delta/24,Qprism.to_numpy(), color='tab:gray', linewidth=2)
    #     axs[1,0].set_ylabel('$Q_{prism}$ (5km)\n$[10^{3}\ m^{3}s^{-1}]$')
    #     axs[1,0].set_yticks(ticks=[20,50,80])
    #     axs[1,1].set_yticks(ticks=[20,50,80])
    #     axs[1,2].set_yticks(ticks=[20,50,80])
    #     axs[1,0].set_ylim(20,80)
    #     axs[1,1].set_ylim(20,80)
    #     axs[1,2].set_ylim(20,80)
    #     # ax0.set_xlim(pd.Timestamp('2020-09-01'), pd.Timestamp('2020-12-31'))
    #     snmid=(np.max(Qprism)+np.min(Qprism))/2
    #     snbg=np.where(Qprism.to_numpy()>snmid, 1, 0)
    #     axs[0,0].pcolor(ot_hours_delta/24, axs[0,0].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True) #slight change to the shading
    #     axs[0,1].pcolor(ot_hours_delta/24, axs[0,1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
    #     axs[0,2].pcolor(ot_hours_delta/24, axs[0,1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
    #     axs[1,0].pcolor(ot_hours_delta/24, axs[1,0].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
    #     axs[1,1].pcolor(ot_hours_delta/24, axs[1,1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
    #     axs[1,2].pcolor(ot_hours_delta/24, axs[0,1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
    #     axs[1,0].grid(True)
    #     axs[1,1].grid(True)
    #     axs[1,2].grid(True)

    # axs[0,0].plot(time_hours/24, zfun.lowpass(count_est, f='godin'), color=linecolor, label=silllenlabel+' total')
    # axs[0,0].plot(time_hours/24, zfun.lowpass(count_est_strict, f='godin'), '--', color=linecolor, label=silllenlabel+' no return')
    # axs[0,1].plot(time_hours/24, zfun.lowpass(count_insill, f='godin'), color=linecolor, label=silllenlabel+' total') 
    # axs[0,1].plot(time_hours/24, zfun.lowpass(count_insill_strict, f='godin'), '--', color=linecolor, label=silllenlabel+' no return') 
    # axs[0,2].plot(time_hours/24, zfun.lowpass(count_in, f='godin'), color=linecolor, label=silllenlabel+' total') 
    # axs[0,2].plot(time_hours/24, zfun.lowpass(count_in_strict, f='godin'), '--', color=linecolor, label=silllenlabel+' no return') 

    axs[i].stackplot(time_hours/24,ret_count_0_ta,ret_count_1_ta,ret_count_2_ta,ret_count_3_ta,ret_count_4plus_ta,colors=['tab:grey',plt.cm.tab20b(2),'tab:cyan','tab:olive','tab:pink'],labels=['0','1','2','3','4+'])
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
# axs[0,0].set_ylabel('Particles remaining in inner basin') #TRY WITH TOTAL PARTICLE COUNT
# axs[0,0].set_title('Released in inner basin below sill height')
# axs[0,0].grid(True)
# axs[0,0].set_xlim(0,120)
# # axs[0,0].set_ylim(0,100)
# axs[0,0].set_ylim(0,count_start_inlow_stay_in[0])

# # axs[0,1].set_xlabel('Days')
# # axs[0,1].set_ylabel('% of particles')
# axs[0,1].set_title('Released in inner basin above sill height')
# axs[0,1].grid(True)
# axs[0,1].set_xlim(0,120)
# # axs[0,1].set_ylim(0,100)
# axs[0,1].set_ylim(0,count_start_inup_stay_in[0])

# axs[1,0].set_xlabel('Days')
# # axs[1,0].set_ylabel('% of particles remaining in outer basin')
# # axs[1,0].set_ylabel('Particles remaining in outer basin') #TRY WITH TOTAL PARTICLE COUNT
# axs[1,0].set_ylabel('Particles remaining') #TRY WITH TOTAL PARTICLE COUNT
# axs[1,0].set_title('Released in outer basin below sill height')
# axs[1,0].grid(True)
# axs[1,0].set_xlim(0,120)
# # axs[1,0].set_ylim(0,100)
# axs[1,0].set_ylim(0,count_start_outlow_stay_out[0])

# axs[1,1].set_xlabel('Days')
# # axs[1,1].set_ylabel('% of particles')
# axs[1,1].set_title('Released in outer basin above sill height')
# axs[1,1].grid(True)
# axs[1,1].set_xlim(0,120)
# # axs[1,1].set_ylim(0,100)
# axs[1,1].set_ylim(0,count_start_outup_stay_out[0])
# axs[1,1].legend(loc='upper right')

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

axs[0].grid(True)
axs[1].grid(True)
axs[2].grid(True)
axs[3].grid(True)
axs[4].grid(True)
# axs[0,].grid(True)
# axs[0,2].grid(True)

axs[0].set_xlabel('Time [days]')
axs[1].set_xlabel('Time [days]')
axs[2].set_xlabel('Time [days]')
axs[3].set_xlabel('Time [days]')
axs[4].set_xlabel('Time [days]')
# axs[1,1].set_xlabel('Days')
# axs[1,2].set_xlabel('Days')

axs[0].set_xlim(0,120)
axs[1].set_xlim(0,120)
axs[2].set_xlim(0,120)
axs[3].set_xlim(0,120)
axs[4].set_xlim(0,120)
# axs[0,1].set_xlim(0,120)
# axs[0,2].set_xlim(0,120)
# axs[1,0].set_xlim(0,120)
# axs[1,1].set_xlim(0,120)
# axs[1,2].set_xlim(0,120)

axs[0].set_ylim(0,18000)
# axs[1].set_ylim(0,18000)
# axs[2].set_ylim(0,18000)
# axs[3].set_ylim(0,18000)
# axs[4].set_ylim(0,18000)
# axs[0,1].set_ylim(0,21000)
# axs[0,2].set_ylim(0,21000)

# axs[0,0].set_title('Whole estuary')
# axs[0,1].set_title('Inner basin + sill')
# axs[0,2].set_title('Inner basin')

axs[0].set_title('5 km',color='tab:red')
axs[1].set_title('10 km',color='tab:orange')
axs[2].set_title('20 km',color='tab:green')
axs[3].set_title('40 km',color='tab:blue')
axs[4].set_title('80 km',color='tab:purple')

axs[0].set_ylabel('Particles remaining')

axs[0].legend(title='Re-entries')
axs[1].legend(title='Re-entries')
axs[2].legend(title='Re-entries')
axs[3].legend(title='Re-entries')
axs[4].legend(title='Re-entries')
# fig.legend(ncol=1, loc='outside right center')

axs[0].set_box_aspect(1)
axs[1].set_box_aspect(1)
axs[2].set_box_aspect(1)
axs[3].set_box_aspect(1)
axs[4].set_box_aspect(1)

axs[0].text(.02, .98, 'A', horizontalalignment='left', verticalalignment='top', transform=axs[0].transAxes, fontsize=14, fontweight='bold')
axs[1].text(.02, .98, 'B', horizontalalignment='left', verticalalignment='top', transform=axs[1].transAxes, fontsize=14, fontweight='bold')
axs[2].text(.02, .98, 'C', horizontalalignment='left', verticalalignment='top', transform=axs[2].transAxes, fontsize=14, fontweight='bold')
axs[3].text(.02, .98, 'D', horizontalalignment='left', verticalalignment='top', transform=axs[3].transAxes, fontsize=14, fontweight='bold')
axs[4].text(.02, .98, 'E', horizontalalignment='left', verticalalignment='top', transform=axs[4].transAxes, fontsize=14, fontweight='bold')


# axs[0,2].legend(loc='upper right')
#legend
# handles, labels = axs[0,0].get_legend_handles_labels()
# handles_reorder = np.concatenate((handles[::2],handles[1::2]),axis=0)
# labels_reorder = np.concatenate((labels[::2],labels[1::2]),axis=0)
# axs[0,0].legend(handles_reorder,labels_reorder,loc='upper right',ncol=2)


fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_tracker2_reentry_count.png'
plt.savefig(fn_fig)
plt.close()
#plt.show()

