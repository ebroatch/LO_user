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
# fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(20,6))
fig, axs = plt.subplots(2,2,figsize=(15,8), gridspec_kw={'height_ratios': [6,1]}) #add sharex later

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

    if i==0:
        sect_name='b3'
        pad=36
        tef_df, vn_list, vec_list = tef_fun.get_two_layer(tef_5km_fn, sect_name)
        tef_df['Q_prism']=tef_df['qprism']/1000
        Qprism = tef_df['Q_prism'].loc['2020-09-04':'2020-12-28'] #cut off extra pad because qprism uses two godin filters
        ot=tef_df.loc['2020-09-04':'2020-12-28'].index
        ot_hours_delta = (((ot - datetime.datetime(2020,9,1,0,0,0)).total_seconds())/3600).to_numpy()
        axs[0,0].set_ylim(0,100)
        axs[0,1].set_ylim(0,100)
        axs[1,0].plot(ot_hours_delta/24,Qprism.to_numpy(), color='tab:gray', linewidth=2) #cut off the weird ends
        axs[1,1].plot(ot_hours_delta/24,Qprism.to_numpy(), color='tab:gray', linewidth=2)
        axs[1,0].set_ylabel('$Q_{prism}$ (5km)\n$[10^{3}\ m^{3}s^{-1}]$')
        axs[1,0].set_yticks(ticks=[20,50,80])
        axs[1,1].set_yticks(ticks=[20,50,80])
        axs[1,0].set_ylim(20,80)
        axs[1,1].set_ylim(20,80)
        # ax0.set_xlim(pd.Timestamp('2020-09-01'), pd.Timestamp('2020-12-31'))
        snmid=(np.max(Qprism)+np.min(Qprism))/2
        snbg=np.where(Qprism.to_numpy()>snmid, 1, 0)
        axs[0,0].pcolor(ot_hours_delta/24, axs[0,0].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True) #slight change to the shading
        axs[0,1].pcolor(ot_hours_delta/24, axs[0,1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
        axs[1,0].pcolor(ot_hours_delta/24, axs[1,0].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
        axs[1,1].pcolor(ot_hours_delta/24, axs[1,1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
        axs[1,0].grid(True)
        axs[1,1].grid(True)

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


    lon_vals = dsr.lon.values
    time_hours = dsr.Time.values
    dsr.close()
    # print('got lon_vals and time\n')
    lon_in = lon_vals >= sillland #these are particles in the inner basin at any given point in time
    lon_in_start_in = lon_in * lon_in[np.newaxis, 0, :]
    # print('got lon_in\n')
    lon_out = (lon_vals <= sillsea) & (lon_vals >= 0)
    lon_out_start_out = lon_out * lon_out[np.newaxis, 0, :]
    # print('got lon_out\n')


    par_in = np.sum(lon_in_start_in,axis=1)
    par_out = np.sum(lon_out_start_out,axis=1)
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
    


    axs[0,0].plot(time_hours/24, zfun.lowpass((par_out/par_out[0])*100, f='godin'), color=linecolor, label=silllenlabel+' tidal avg') #PLOT ONLY OUTER PARTICLES REMAINING
    axs[0,1].plot(time_hours/24, zfun.lowpass((par_in/par_in[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY INNER PARTICLES REMAINING    
    # ax1.plot(time_hours/24, (par_out/par_out[0])*100, color=linecolor, label=silllenlabel) #TRY WITH NO FILTERING
    # ax3.plot(time_hours/24, (par_in/par_in[0])*100, color=linecolor, label=silllenlabel)
    # print('plotted\n')
    # sys.stdout.flush()
    
    #dsg.close()

    #ADD CALCULATIONS TO PRINT FOR EACH MODEL
    pad=36 #could use 35 but 36 is nice bc it gives exactly 1.5 days removed on each end
    
    par_out_ta = zfun.lowpass(par_out, f='godin')[pad:-pad+1]
    par_in_ta = zfun.lowpass(par_in, f='godin')[pad:-pad+1]

    par_out_frac_ta = zfun.lowpass((par_out/par_out[0]), f='godin')[pad:-pad+1]
    par_in_frac_ta = zfun.lowpass((par_in/par_in[0]), f='godin')[pad:-pad+1]

    t_ta = time_hours[pad:-pad+1]
    par_scale_out = par_out_ta[0]
    par_scale_in = par_in_ta[0]
    t_scale = t_ta[-1]
    t_frac = t_ta / t_scale

    par_out_frac_raw = par_out/par_out[0]
    par_in_frac_raw = par_in/par_in[0]
    t_scale_raw = time_hours[-1]
    t_frac_raw = time_hours / t_scale_raw

    #quick calc
    T_e_out = -(t_ta[-1] / np.log(par_out_frac_ta[-1]/1))/24 #e-folding time in days from tidally averaged values (negative because it is exponential decrease)
    T_e_in = -(t_ta[-1] / np.log(par_in_frac_ta[-1]/1))/24
    print('\nT_e_out:\n')
    print(T_e_out)
    print('\nT_e_in:\n')
    print(T_e_in)
    T_e_out_raw = -(time_hours[-1] / np.log(par_out[-1]/par_out[0]))/24 #e-folding time in days from raw number of particles
    T_e_in_raw = -(time_hours[-1] / np.log(par_in[-1]/par_in[0]))/24
    # print('quick e-folding calc done\n')
    sys.stdout.flush()

    #exponential fit with offset
    def func(x, a, b, c):
        return a * np.exp(-b * x) + c

    p0=(1,1,0)
    popt_out, pcov_out = curve_fit(func, t_frac, par_out_frac_ta, p0=p0)
    popt_in, pcov_in = curve_fit(func, t_frac, par_in_frac_ta, p0=p0)

    T_e_out_fit = (1/popt_out[1])*t_scale/24 #T_e from fit in days
    T_e_in_fit = (1/popt_in[1])*t_scale/24 #T_e from fit in days   

    #exponential fit
    def func2(x, a, b):
        return a * np.exp(-x/b) #change b to be time units

    p02=(1,1)
    popt_out2, pcov_out2 = curve_fit(func2, t_frac, par_out_frac_ta, p0=p02)
    popt_in2, pcov_in2 = curve_fit(func2, t_frac, par_in_frac_ta, p0=p02)

    T_e_out_fit2 = (popt_out2[1])*t_scale/24 #T_e from two parameter fit in days
    T_e_in_fit2 = (popt_in2[1])*t_scale/24 #T_e from two parameter fit in days  
    A_out_fit2 = popt_out2[0]
    A_in_fit2 = popt_in2[0] 

    print('\nT_e_out from two-param fit:\n')
    print(T_e_out_fit2)
    print('\nT_e_in from two-param fit:\n')
    print(T_e_in_fit2)

    #exponential fit of raw data
    p03=(1,1)
    popt_out3, pcov_out3 = curve_fit(func2, t_frac_raw, par_out_frac_raw, p0=p03)
    popt_in3, pcov_in3 = curve_fit(func2, t_frac_raw, par_in_frac_raw, p0=p03)

    T_e_out_fit3 = (popt_out3[1])*t_scale_raw/24 #T_e from two parameter fit in days
    T_e_in_fit3 = (popt_in3[1])*t_scale_raw/24 #T_e from two parameter fit in days
    A_out_fit3 = popt_out3[0]
    A_in_fit3 = popt_in3[0]

    print('\nT_e_out from two-param fit (unfiltered data):\n')
    print(T_e_out_fit3)
    print('\nT_e_in from two-param fit (unfiltered data):\n')
    print(T_e_in_fit3)

    #add fits to plot
    # axs[0,0].plot(time_hours/24, A_out_fit3 * np.exp((-time_hours/24)/T_e_out_fit3) * 100, linestyle = '--', color=linecolor2, label=silllenlabel+' fit')
    # axs[0,1].plot(time_hours/24, A_in_fit3 * np.exp((-time_hours/24)/T_e_in_fit3) * 100, linestyle = '--', color=linecolor2, label=silllenlabel)

    axs[0,0].plot(time_hours/24, A_out_fit2 * np.exp((-time_hours/24)/T_e_out_fit2) * 100, linestyle = '--', color=linecolor2, label=silllenlabel+' fit')
    axs[0,1].plot(time_hours/24, A_in_fit2 * np.exp((-time_hours/24)/T_e_in_fit2) * 100, linestyle = '--', color=linecolor2, label=silllenlabel)


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

axs[0,0].set_ylabel('% of particles')
axs[0,0].set_title('Particles released in outer basin')
# ax1.set_title('Particles released in outer basin (unfiltered)')
#ax1.legend(loc='best')
axs[0,0].grid(True)
axs[0,0].set_xlim(0,120)
axs[0,0].set_ylim(0,100)

#legend
handles, labels = axs[0,0].get_legend_handles_labels()
handles_reorder = np.concatenate((handles[::2],handles[1::2]),axis=0)
labels_reorder = np.concatenate((labels[::2],labels[1::2]),axis=0)
axs[0,0].legend(handles_reorder,labels_reorder,loc='upper right',ncol=2)

# ax2.set_xlabel('Days')
# #ax1.set_ylabel('Number of particles')
# ax2.set_title('Particles released on sill')
# #ax2.legend(loc='best')
# ax2.grid(True)
# ax2.set_xlim(0,120)
# ax2.set_ylim(0,100)

# axs[0,1].set_xlabel('Days')
#ax3.set_ylabel('Number of particles')
axs[0,1].set_title('Particles released in inner basin')
# ax3.set_title('Particles released in inner basin (unfiltered)')
#ax3.legend(loc='best')
axs[0,1].grid(True)
axs[0,1].set_xlim(0,120)
axs[0,1].set_ylim(0,100)

axs[1,0].set_xlim(0,120)
axs[1,1].set_xlim(0,120)
axs[1,0].set_ylim(20,80)
axs[1,1].set_ylim(20,80)
axs[1,0].set_xlabel('Days')
axs[1,1].set_xlabel('Days')

fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_tracker2_fitandcalc_sn.png'
plt.savefig(fn_fig)
plt.close()
#plt.show()
