"""
Plot results of a particle tracking experiment.
"""


import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from scipy.optimize import curve_fit

from lo_tools import Lfun
from lo_tools import zfun
from lo_tools import plotting_functions as pfun
Ldir = Lfun.Lstart()

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

# # subsample output for plotting #SKIP SUBSAMPLING
# npmax = 600 # max number of points to plot
# step = max(1,int(np.floor(NP/npmax)))

# lon = dsr.lon.values[:,::step]
# lat = dsr.lat.values[:,::step]
lon = dsr.lon
lonest = (lon>0)
lonest = lonest.astype(int)
partest = lonest.sum(dim='Particle')
partest_ta = zfun.lowpass(partest.values,f='godin',nanpad=True)[35:-35]
tplot = partest.Time.values[35:-35]

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

p0=(30000,-0.0005,1)
#popt, pcov = curve_fit(func, tplot, partest_ta, p0=p0)
popt, pcov = curve_fit(func, partest.Time, partest, p0=p0)
#pfit=np.polyfit(tplot,np.log(partest_ta),1)

# PLOTTING - SPAGHETTI PLOT
plt.close('all')
fig, ax = plt.subplots(1,1)

ax.set_xlabel('Hours')
ax.set_ylabel('Particles in estuary')
ax.set_title('Whole estuary residence time')
# add the tracks (packed [time, particle])
# regular spaghetti plots
ax.plot(partest.Time, partest, '.c', label='Raw particle #') #data
ax.plot(tplot, partest_ta, '-b', label='Tidally averaged') #data
# plt.plot(tplot, func(tplot, *popt), '--r', label='fit: a=%5.3f, b=%5.3f' % tuple(popt)) #fit
# plt.plot(tplot, func(tplot, *popt), '--r', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)) #fit
plt.plot(partest.Time, func(partest.Time, *popt), '--r', label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt)) #fit
# plt.plot(tplot, np.exp(pfit[1])*np.exp(pfit[0]*tplot), '--r', label='Fit')
ax.legend(loc='best')


#plt.show()
pfun.end_plot()

# #PLOTTING - HISTOGRAMS
# fig, axs = plt.subplots(5,1,sharex=True)
# for j in range(5):
#     hour=j*180
#     axs[j].set_title('t='+str(hour)+'h')
#     axs[j].hist(dsr['lon'].sel(Time=hour),bins=20,alpha=0.5)
#     #axs[j].set_ylim(0, 30)

#plt.show()
fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rt.png'
plt.savefig(fn_fig)

dsr.close()
dsg.close()