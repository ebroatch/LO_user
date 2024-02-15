"""
Plot results of a particle tracking experiment.
"""
from lo_tools import Lfun
from lo_tools import zfun
from lo_tools import plotting_functions as pfun
from lo_user_tools import llxyfun
Ldir = Lfun.Lstart()

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from scipy.optimize import curve_fit

plt.close('all')
fig, [ax1,ax2,ax3] = plt.subplots(1,3,figsize=(20,6))

for i in range(3):
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

    # # subsample output for plotting #SKIP SUBSAMPLING
    # npmax = 600 # max number of points to plot
    # step = max(1,int(np.floor(NP/npmax)))

    # lon = dsr.lon.values[:,::step]
    # lat = dsr.lat.values[:,::step]
    # lon = dsr.lon

    if i==0:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(45e3,0,45)
        linst = ':'
    elif i==1:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(60e3,0,45)
        linst = '-'
    elif i==2:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(120e3,0,45)
        linst = '--'


    lon1 = dsr.lon.where((dsr.lon.sel(Time=0)<sillsea),drop=True)
    lon2 = dsr.lon.where((dsr.lon.sel(Time=0)>=sillsea) & (dsr.lon.sel(Time=0)<sillland),drop=True)
    lon3 = dsr.lon.where((dsr.lon.sel(Time=0)>=sillland),drop=True)

    par1_ocn=(lon1<0).astype(int).sum(dim='Particle')
    par1_out=((lon1>=0) & (lon1<sillsea)).astype(int).sum(dim='Particle')
    par1_sill=((lon1>=sillsea) & (lon1<sillland)).astype(int).sum(dim='Particle')
    par1_in=(lon1>=sillland).astype(int).sum(dim='Particle')

    par2_ocn=(lon2<0).astype(int).sum(dim='Particle')
    par2_out=((lon2>=0) & (lon2<sillsea)).astype(int).sum(dim='Particle')
    par2_sill=((lon2>=sillsea) & (lon2<sillland)).astype(int).sum(dim='Particle')
    par2_in=(lon2>=sillland).astype(int).sum(dim='Particle')

    par3_ocn=(lon3<0).astype(int).sum(dim='Particle')
    par3_out=((lon3>=0) & (lon3<sillsea)).astype(int).sum(dim='Particle')
    par3_sill=((lon3>=sillsea) & (lon3<sillland)).astype(int).sum(dim='Particle')
    par3_in=(lon3>=sillland).astype(int).sum(dim='Particle')

    # lonest = (lon>0)
    # lonest = lonest.astype(int)
    # partest = lonest.sum(dim='Particle')
    # partest_ta = zfun.lowpass(partest.values,f='godin',nanpad=True)[35:-35]
    # tplot = partest.Time.values[35:-35]


    ax1.plot(par1_ocn.Time/24, zfun.lowpass((par1_ocn/par1_out.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:green', label='Ocean')
    ax1.plot(par1_out.Time/24, zfun.lowpass((par1_out/par1_out.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:cyan', label='Outer basin')
    ax1.plot(par1_sill.Time/24, zfun.lowpass((par1_sill/par1_out.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:purple', label='Sill')
    ax1.plot(par1_in.Time/24, zfun.lowpass((par1_in/par1_out.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:pink', label='Inner basin')

    #ax.set_ylim(20000,35000)


    ax2.plot(par2_ocn.Time/24, zfun.lowpass((par2_ocn/par2_sill.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:green', label='Ocean')
    ax2.plot(par2_out.Time/24, zfun.lowpass((par2_out/par2_sill.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:cyan', label='Outer basin')
    ax2.plot(par2_sill.Time/24, zfun.lowpass((par2_sill/par2_sill.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:purple', label='Sill')
    ax2.plot(par2_in.Time/24, zfun.lowpass((par2_in/par2_sill.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:pink', label='Inner basin')


    ax3.plot(par3_ocn.Time/24, zfun.lowpass((par3_ocn/par3_in.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:green', label='Ocean')
    ax3.plot(par3_out.Time/24, zfun.lowpass((par3_out/par3_in.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:cyan', label='Outer basin')
    ax3.plot(par3_sill.Time/24, zfun.lowpass((par3_sill/par3_in.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:purple', label='Sill')
    ax3.plot(par3_in.Time/24, zfun.lowpass((par3_in/par3_in.sel(Time=0)).values, f='godin'), linestyle=linst, color='tab:pink', label='Inner basin')

    dsr.close()
    dsg.close()


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
ax1.set_ylabel('%% of particles')
ax1.set_title('Particles released in outer basin')
ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlim(0,120)
ax1.set_ylim(0,100)

ax2.set_xlabel('Days')
#ax1.set_ylabel('Number of particles')
ax2.set_title('Particles released on sill')
ax2.legend(loc='best')
ax2.grid(True)
ax2.set_xlim(0,120)
ax2.set_ylim(0,100)

ax3.set_xlabel('Days')
#ax3.set_ylabel('Number of particles')
ax3.set_title('Particles released in inner basin')
ax3.legend(loc='best')
ax3.grid(True)
ax3.set_xlim(0,120)
ax3.set_ylim(0,100)


fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_osm.png'
plt.savefig(fn_fig)
plt.close()
#plt.show()

