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

plt.close('all')
# fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(20,6))
fig, [ax1,ax3] = plt.subplots(1,2,figsize=(15,6))

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


    lon_vals = dsr.lon.values
    time_hours = dsr.Time.values
    dsr.close()
    print('got lon_vals and time\n')
    lon_in = lon_vals >= sillland #these are particles in the inner basin at any given point in time
    lon_in_start_in = lon_in * lon_in[np.newaxis, 0, :]
    print('got lon_in\n')
    lon_out = (lon_vals <= sillsea) & (lon_vals >= 0)
    lon_out_start_out = lon_out * lon_out[np.newaxis, 0, :]
    print('got lon_out\n')


    par_in = np.sum(lon_in_start_in,axis=1)
    par_out = np.sum(lon_out_start_out,axis=1)
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
    


    ax1.plot(time_hours/24, zfun.lowpass((par_out/par_out[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY OUTER PARTICLES REMAINING
    ax3.plot(time_hours/24, zfun.lowpass((par_in/par_in[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY INNER PARTICLES REMAINING    
    # ax1.plot(time_hours/24, (par_out/par_out[0])*100, color=linecolor, label=silllenlabel) #TRY WITH NO FILTERING
    # ax3.plot(time_hours/24, (par_in/par_in[0])*100, color=linecolor, label=silllenlabel)
    print('plotted\n')
    sys.stdout.flush()
    
    #dsg.close()

    #ADD CALCULATIONS TO PRINT FOR EACH MODEL
    pad=36 #could use 35 but 36 is nice bc it gives exactly 1.5 days removed on each end
    par_out_frac_ta = zfun.lowpass((par_out/par_out[0]), f='godin')[pad:-pad+1]
    par_in_frac_ta = zfun.lowpass((par_in/par_in[0]), f='godin')[pad:-pad+1]

    par_out_ta = zfun.lowpass(par_out, f='godin')[pad:-pad+1]
    par_in_ta = zfun.lowpass(par_in, f='godin')[pad:-pad+1]
    t_ta = time_hours[pad:-pad+1]

    #quick calc
    t_e_out = (t_ta[-1] / np.log(par_out_frac_ta[-1]/1))/24 #e-folding time in days from tidally averaged values
    t_e_in = (t_ta[-1] / np.log(par_in_frac_ta[-1]/1))/24

    t_e_out_raw = (time_hours[-1] / np.log(par_out[-1]/par_out[0]))/24 #e-folding time in days from raw number of particles
    t_e_in_raw = (time_hours[-1] / np.log(par_in[-1]/par_in[0]))/24
    print('quick e-folding calc done\n')
    sys.stdout.flush()

    #exponential fit
    def func(x, a, b, c):
        return a * np.exp(-b * x) + c

    p0=(1,0.0005,0)
    popt_out, pcov_out = curve_fit(func, t_ta, par_out_frac_ta, p0=p0)
    popt_in, pcov_in = curve_fit(func, t_ta, par_in_frac_ta, p0=p0)

    def func2(x, a, b):
        return a * np.exp(-b * x)

    p02=(1,0.0005)
    popt_out2, pcov_out2 = curve_fit(func2, t_ta, par_out_frac_ta, p0=p02)
    popt_in2, pcov_in2 = curve_fit(func2, t_ta, par_in_frac_ta, p0=p02)

    #fitting




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
ax1.set_xlabel('Days')
ax1.set_ylabel('% of particles')
ax1.set_title('Particles released in outer basin')
# ax1.set_title('Particles released in outer basin (unfiltered)')
#ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlim(0,120)
ax1.set_ylim(0,100)

# ax2.set_xlabel('Days')
# #ax1.set_ylabel('Number of particles')
# ax2.set_title('Particles released on sill')
# #ax2.legend(loc='best')
# ax2.grid(True)
# ax2.set_xlim(0,120)
# ax2.set_ylim(0,100)

ax3.set_xlabel('Days')
#ax3.set_ylabel('Number of particles')
ax3.set_title('Particles released in inner basin')
# ax3.set_title('Particles released in inner basin (unfiltered)')
#ax3.legend(loc='best')
ax3.grid(True)
ax3.set_xlim(0,120)
ax3.set_ylim(0,100)
ax3.legend(loc='upper right')


fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_tracker2_fitandcalc.png'
plt.savefig(fn_fig)
plt.close()
#plt.show()
