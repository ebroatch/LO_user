"""
Plot average salinity of estuary and inner basin over year run and four months used for analysis

To test on mac:
run bulk_plot -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31 -test True


"""
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pickle
from time import time
import pandas as pd
import xarray as xr

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import tef_fun
import scipy.signal

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

silllens_plot=[5,10,20,40,80]

sect_est ='a1'
sect_in='b5'

silllens=['5km', '10km', '20km', '40km', '80km']
gridnames = ['sill5km', 'sill10km', 'sill20kmdeep', 'sill40km', 'sill80km']
gctags=['sill5km_c0', 'sill10km_c0', 'sill20kmdeep_c0', 'sill40km_c0', 'sill80km_c0']
gtagexs=['sill5km_t0_xa0', 'sill10km_t2_xa0', 'sill20kmdeep_t2_xa0', 'sill40km_t2_xa0', 'sill80km_t2_xa0']
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
plot_colors = ['tab:red','tab:orange','tab:green','tab:blue','tab:purple']

fig, axs = plt.subplots(2,2,figsize=(15,8))

#Loop over sill lengths
for i in range(len(gctags)):
# for i in range(len(gctags)-1):
# for i in [0,2]:
    #model and extraction info
    print(silllens[i])
    gctag=gctags[i]
    gtagex=gtagexs[i]
    gridname=gridnames[i]
    in_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('bulk_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
    seg_est_inner_fn = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('segments_' + Ldir['ds0'] + '_' + Ldir['ds1'] + '_' + gridname + '_cei_rivA1.nc')
    
    # Next, get the average salinity of the estuary and inner basin
    # These use the segment extractions for collections cei
    # the inner basin is segment b5i_p
    # for the whole estuary, need to get the weighted average of segment b5i_p and b5i_m
    seg_est_inner_ds = xr.open_dataset(seg_est_inner_fn)
    sbar_inner_ts = seg_est_inner_ds.salt
    sbar_outer_sill_ts = seg_est_inner_ds.salt
    vol_inner_ts = seg_est_inner_ds.volume
    vol_outer_sill_ts = seg_est_inner_ds.volume
    sbar_est_ts = ((sbar_inner_ts*vol_inner_ts)+(sbar_outer_sill_ts*vol_outer_sill_ts))/(vol_inner_ts+vol_outer_sill_ts)
    # get the time average values
    sbar_inner = sbar_inner_ts.mean()
    sbar_est = sbar_est_ts.mean()

    #plot the salinities #MIGHT NEED TO TIDALLY AVERAGE??
    axs[0,0].plot(sbar_est_ts,c=plot_colors[i],label=silllens[i])
    axs[0,1].plot(sbar_inner_ts,c=plot_colors[i],label=silllens[i])
    axs[1,0].plot(sbar_est_ts,c=plot_colors[i],label=silllens[i]) #we will change the axis limits to make these more zoomed
    axs[1,1].plot(sbar_inner_ts,c=plot_colors[i],label=silllens[i])

axs[0,0].set_title(r'Whole estuary mean salinity $\bar{s}$')
axs[1,0].set_title(r'Inner basin mean salinity $\bar{s}$')
axs[0,0].text(.04, .02, 'A', horizontalalignment='left', verticalalignment='bottom', transform=axs[0,0].transAxes, fontsize=14, fontweight='bold')
axs[0,1].text(.04, .02, 'B', horizontalalignment='left', verticalalignment='bottom', transform=axs[0,1].transAxes, fontsize=14, fontweight='bold')
axs[1,0].text(.04, .02, 'C', horizontalalignment='left', verticalalignment='bottom', transform=axs[1,0].transAxes, fontsize=14, fontweight='bold')
axs[1,1].text(.04, .02, 'D', horizontalalignment='left', verticalalignment='bottom', transform=axs[1,1].transAxes, fontsize=14, fontweight='bold')
axs[0,0].xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
axs[0,1].xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
axs[1,0].xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
axs[1,1].xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
axs[0,0].set_xlim('2020-01-01','2020-12-31')
axs[0,1].set_xlim('2020-01-01','2020-12-31')
axs[1,0].set_xlim('2020-09-01','2020-12-31')
axs[1,1].set_xlim('2020-09-01','2020-12-31')
axs[0,0].set_ylim(0,34)
axs[0,1].set_ylim(0,34)
axs[1,0].set_ylim(25,34)
axs[1,1].set_ylim(25,34)
axs[0,0].grid(True)
axs[0,1].grid(True)
axs[1,0].grid(True)
axs[1,1].grid(True)
axs[0,0].set_ylabel('Salinity [psu]')
axs[1,0].set_ylabel('Salinity [psu]')
axs[1,0].set_xlabel('Date')
axs[1,1].set_xlabel('Date')
axs[0,0].legend(loc='upper left')
axs[0,1].legend(loc='upper left')
axs[1,0].legend(loc='upper left')
axs[1,1].legend(loc='upper left')

fn_fig = Ldir['LOo'] / 'plots' / 'salinity_spinup_estuary_inner.png'
plt.savefig(fn_fig)
plt.close()

# #plot
# ax1.plot(silllens_plot,flushing_est_days,c='k',marker='o',label='Volume flushing')
# ax1.plot(silllens_plot,flushing_fresh_est_days,c='tab:cyan',marker='o',label='Freshwater flushing')
# ax1.plot(silllens_plot,flushing_salt_est_days,c='tab:olive',marker='o',label='Saltwater flushing')
# ax2.plot(silllens_plot,flushing_inner_days,c='k',marker='o',label='Volume flushing')
# ax2.plot(silllens_plot,flushing_fresh_inner_days,c='tab:cyan',marker='o',label='Freshwater flushing')
# ax2.plot(silllens_plot,flushing_salt_inner_days,c='tab:olive',marker='o',label='Saltwater flushing')

# #add plot elements
# ax1.set_xlabel('Sill length')
# ax1.set_ylabel('Flushing time [days]')
# ax1.set_ylim(0,85)
# ax1.set_ylim(30,75)
# ax1.set_title('Estuary')
# ax1.grid(True)
# ax1.legend()

# ax2.set_xlabel('Sill length')
# # ax2.set_ylabel('Flushing time [days]')
# ax2.set_ylim(0,85)
# ax2.set_ylim(30,75)
# ax2.set_title('Inner basin')
# ax2.grid(True)
# # ax2.legend()
# plt.suptitle('Volume, freshwater, and saltwater flushing times')

# # h, l = ax.get_legend_handles_labels()
# # ph = [plt.plot([],marker="", ls="")[0]]*2
# # handles = ph + h
# # labels = [r'Inner basin reflux:', r'Outer basin reflux:'] + l
# # # ax.legend(handles, labels, ncol=6)
# # ax.legend(handles,labels,ncol=5)

# fn_fig = Ldir['LOo'] / 'plots' / 'flushing_times.png' #UNCOMMENT TO PLOT
# plt.savefig(fn_fig)
# plt.close()

#update plot style for thesis
# pfun.start_plot(fs=14)
# fig, [ax1,ax2] = plt.subplots(1,2,figsize=(12,6))
# fig, [ax1,ax2] = plt.subplots(1,2,figsize=(8,5))
#Only plot the reflux fractions
# ax1.plot(silllens_plot,flushing_est_days,c='k',marker='o',label='Volume')
# ax1.plot(silllens_plot,flushing_fresh_est_days,c='tab:cyan',marker='o',label='Freshwater')
# ax1.plot(silllens_plot,flushing_salt_est_days,c='tab:olive',marker='o',label='Saltwater')
# ax2.plot(silllens_plot,flushing_inner_days,c='k',marker='o',label='Volume')
# ax2.plot(silllens_plot,flushing_fresh_inner_days,c='tab:cyan',marker='o',label='Freshwater')
# ax2.plot(silllens_plot,flushing_salt_inner_days,c='tab:olive',marker='o',label='Saltwater')
# ax1.plot(silllens_plot,flushing_est_days,c='teal',marker='o',label='Volume')
# ax1.plot(silllens_plot,flushing_fresh_est_days,c='lightskyblue',marker='o',label='Freshwater')
# ax1.plot(silllens_plot,flushing_salt_est_days,c='tab:olive',marker='o',label='Saltwater')
# ax2.plot(silllens_plot,flushing_inner_days,c='teal',marker='o',label='Volume')
# ax2.plot(silllens_plot,flushing_fresh_inner_days,c='lightskyblue',marker='o',label='Freshwater')
# ax2.plot(silllens_plot,flushing_salt_inner_days,c='tab:olive',marker='o',label='Saltwater')
# ax1.set_xlabel('Sill length [km]')
# ax2.set_xlabel('Sill length [km]')
# ax1.set_ylabel('Flushing time [days]')
# ax1.set_title('Whole estuary')
# ax2.set_title('Inner basin')
# ax1.grid(True)
# ax2.grid(True)
# ax1.legend(loc='lower right')
# ax2.legend(loc='lower right')
# ax1.set_xlim(0,85)
# ax2.set_xlim(0,85)
# ax1.set_ylim(30,75)
# ax2.set_ylim(30,75)
# ax1.set_box_aspect(1)
# ax2.set_box_aspect(1)
# ax1.text(.02, .98, 'A', horizontalalignment='left', verticalalignment='top', transform=ax1.transAxes, fontsize=14, fontweight='bold')
# ax2.text(.02, .98, 'B', horizontalalignment='left', verticalalignment='top', transform=ax2.transAxes, fontsize=14, fontweight='bold')
fn_fig = Ldir['LOo'] / 'plots' / 'flushing_times.png'
plt.savefig(fn_fig)
plt.close()

# print('\nalpha_21 (outer basin reflux): ')
# print(alpha_21_basic)
# print('\nalpha_31 (efflux from outer basin): ')
# print(alpha_31_basic)
# print('\nalpha_34 (inner basin reflux): ')
# print(alpha_34_basic)
# print('\nalpha_24 (efflux from inner basin): ')
# print(alpha_24_basic)

# print('\nalpha_21+alpha_31 (outer basin fractions): ')
# print(alpha_21_basic+alpha_31_basic)
# print('\nalpha_24+alpha_34 (outer basin fractions): ')
# print(alpha_24_basic+alpha_34_basic)

