"""
Plot bulk fluxes as a time series.

To test on mac:
run bulk_plot -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31 -test True


"""
import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time
import pandas as pd
import xarray as xr

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import flux_fun

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('salt_integral_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('si_plots_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

sect_list = ['a1.nc','b3.nc']
sect_label = ['Estuary', 'Inner basin']
plot_color = ['tab:cyan','tab:pink']
plot_color2 = ['tab:blue','tab:red']

plt.close('all')
fs = 14
pfun.start_plot(fs=fs, figsize=(10,15))

fig, axs = plt.subplots(3, 1, sharex=True,figsize=(15,10),gridspec_kw={'height_ratios': [2,2,2]})

for i in range(len(sect_list)):
    sect_name = sect_list[i]
    fn = in_dir / sect_name.replace('.nc','.p')
    #fn = in_dir / sect_name #change to this once filenames are fixed
    ds = xr.open_dataset(fn)

    #Plot salinity
    axs[0].plot(ds.time, ds.s_bar, color=plot_color[i], linewidth=1, label=sect_label[i] + r'$\bar{s}$')

    #Plot tidally averaged salinity
    axs[0].plot(ds.time, zfun.lowpass(ds.s_bar.values,f='godin'), color=plot_color2[i], linewidth=2, label=sect_label[i] + r'$\langle \bar{s} \rangle$')
    
    #Calculate spring-neap p2p and 15 day trend (best-fit slope over 15 days)
    snrange=[]
    tplot=[]
    slope15d=[]
    for j in range(180, len(ds.time)-180):
        smax=np.nanmax(ds.s_bar.values[j-180:j+181])
        smin=np.nanmin(ds.s_bar.values[j-180:j+181])
        tplot.append(ds.time.values[j])
        snrange.append(smax-smin)
        p = np.polyfit(np.arange(0,361),ds.s_bar.values[j-180,j+181],1)
        slope15d.append(p[0])
    tplot=np.asarray(tplot)
    snrange=np.asarray(snrange)
    slope15d=np.asarray(slope15d)

    #Spinup condition #add after look at plot
    #tspin=tplot[np.argmax(np.where(snrange>(10*360*slope15d)))] #find first time where spring-neap range is 10x greater than 15 day trend

    #Plot spring-neap p2p
    axs[1].plot(tplot,snrange, color=plot_color2[i], linewidth=2, label=sect_label[i])
    #axs[1].axvline(x=tspin, color=plot_color[i], linestyle='--', label=sect_label[i] + ' spinup')
    
    #Plot best fit slope * 15 days (360h)
    axs[2].plot(tplot, slope15d*360, color=plot_color2[i], linewidth=2, label=sect_label[i])
    #axs[2].axvline(x=tspin, color=plot_color[i], linestyle='--')

    #SD = pickle.load(open(in_dir / sect_name, 'rb'))

    # FR = SD['FR']
    # FE = SD['FE']
    # FT = SD['FT']
    # FTL = SD['FTL']
    # FTV = SD['FTV']
    # F = SD['F']
    # ot = SD['ot']
                    
    # labels and colors
    # ylab_dict = {'Q': r'Transport $[10^{3}\ m^{3}s^{-1}]$',
    #             'salt': r'Salinity $[g\ kg^{-1}]$'}
    # ylab_dict = {'Q': r'$Q_{in}\ [10^{3}\ m^{3}s^{-1}]$',
    #             'salt': r'$s_{in}\ [g\ kg^{-1}]$',
    #             'deltas': r'$\Delta s\ [g\ kg^{-1}]$'}
    # p_color = 'r'
    # m_color = 'b'
    lw = 2
        
    # fig = plt.figure()
    
    # ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
    # ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
    # ax3 = plt.subplot2grid((1,3), (0,2)) # map
    
    #ot = bulk['ot'] # (same as tef_df.index)
    
    # axs[0].plot(ot,FR, color='tab:blue', linewidth=lw, label=sect_label[i])
    # axs[2].plot(ot,FE, color='tab:green', linewidth=lw, label=sect_label[i])
    # axs[3].plot(ot,FT, color='tab:red', linewidth=lw, label=sect_label[i])
    # axs[2].plot(ot, FR, color=plot_color[i], linewidth=lw, label=sect_label[i])
    # axs[0].plot(ot, FE, color=plot_color[i], linewidth=lw, label=sect_label[i])
    # axs[1].plot(ot, FT, color=plot_color[i], linewidth=lw, label=sect_label[i])
    # axs[i].plot(ot,FTL, linestyle = '--', color='tab:pink', linewidth=lw, label=r'$F_{TL}$')
    # axs[i].plot(ot,FTV, linestyle = ':', color='tab:orange', linewidth=lw, label=r'$F_{TV}$')
axs[0].grid(True)
axs[1].grid(True) 
axs[2].grid(True)   
axs[0].set_ylabel(r'Salinity $[g\ kg^{-1}]$')
#axs[1].set_ylabel(r'$\langle \bar{s} \rangle [g\ kg^{-1}]$')
axs[1].set_ylabel(r'Spring-neap salinity difference $[g\ kg^{-1}]$') #change this
axs[2].set_ylabel(r'15 day salinity trend $[g\ kg^{-1}]$') #change this
# axs[2].set_ylim(-3.2e4,-2.8e4)
# axs[0].set_ylim(-1e4,5e4)
# axs[1].set_ylim(-2e4,4e4)
axs[0].legend(loc='lower right')
axs[2].set_xlim(pd.Timestamp('2020-01-01'), pd.Timestamp('2020-07-31'))
#plt.suptitle('Standard decomposition')
plt.savefig(out_dir / ('si_plot.png'))
plt.close()
pfun.end_plot()