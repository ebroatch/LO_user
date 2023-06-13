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
sillmid = llxyfun.x2lon(44e3,0,45)
lon1 = dsr.lon.where((dsr.lon.sel(Time=0)<sillmid),drop=True)
lon2 = dsr.lon.where((dsr.lon.sel(Time=0)>sillmid),drop=True)

par1_ocn=(lon1<0).astype(int).sum(dim='Particle')
par1_out=((lon1>0) & (lon1<sillmid)).astype(int).sum(dim='Particle')
par1_in=(lon1>sillmid).astype(int).sum(dim='Particle')

par2_ocn=(lon2<0).astype(int).sum(dim='Particle')
par2_out=((lon2>0) & (lon2<sillmid)).astype(int).sum(dim='Particle')
par2_in=(lon2>sillmid).astype(int).sum(dim='Particle')

# lonest = (lon>0)
# lonest = lonest.astype(int)
# partest = lonest.sum(dim='Particle')
# partest_ta = zfun.lowpass(partest.values,f='godin',nanpad=True)[35:-35]
# tplot = partest.Time.values[35:-35]

plt.close('all')
fig, [ax1,ax2] = plt.subplots(1,2,figsize=(16,6))

ax1.set_xlabel('Days')
ax1.set_ylabel('Number of particles')
ax1.set_title('Particles released in outer basin')
ax1.plot(par1_ocn.Time/24, par1_ocn, '-', color='tab:green', label='Ocean')
ax1.plot(par1_out.Time/24, par1_out, '-', color='tab:cyan', label='Outer basin')
ax1.plot(par1_in.Time/24, par1_in, '-', color='tab:pink', label='Inner basin')
ax1.legend(loc='best')
ax1.grid(True)
ax1.set_xlim(0,120)
ax1.set_ylim(0,18000)
#ax.set_ylim(20000,35000)

ax2.set_xlabel('Days')
#ax1.set_ylabel('Number of particles')
ax2.set_title('Particles released in inner basin')
ax2.plot(par2_ocn.Time/24, par2_ocn, '-', color='tab:green', label='Ocean')
ax2.plot(par2_out.Time/24, par2_out, '-', color='tab:cyan', label='Outer basin')
ax2.plot(par2_in.Time/24, par2_in, '-', color='tab:pink', label='Inner basin')
ax2.legend(loc='best')
ax2.grid(True)
ax2.set_xlim(0,120)
ax2.set_ylim(0,18000)

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
fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins.png'
plt.savefig(fn_fig)
plt.close()
#plt.show()

dsr.close()
dsg.close()