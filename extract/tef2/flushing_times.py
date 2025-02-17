"""
Plot bulk fluxes as a time series.

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

flushing_est=np.zeros(5)
flushing_inner=np.zeros(5)
flushing_est_days=np.zeros(5)
flushing_inner_days=np.zeros(5)

# flushing_fresh_est=np.zeros(5)
# flushing_fresh_inner=np.zeros(5)
# flushing_salt_est=np.zeros(5)
# flushing_salt_inner=np.zeros(5)
flushing_fresh_est_days=np.zeros(5)
flushing_fresh_inner_days=np.zeros(5)
flushing_salt_est_days=np.zeros(5)
flushing_salt_inner_days=np.zeros(5)

silllens_plot=[5,10,20,40,80]

sect_est ='a1'
sect_in='b5'

silllens=['5km', '10km', '20km', '40km', '80km']
gridnames = ['sill5km', 'sill10km', 'sill20kmdeep', 'sill40km', 'sill80km']
gctags=['sill5km_c0', 'sill10km_c0', 'sill20kmdeep_c0', 'sill40km_c0', 'sill80km_c0']
gtagexs=['sill5km_t0_xa0', 'sill10km_t2_xa0', 'sill20kmdeep_t2_xa0', 'sill40km_t2_xa0', 'sill80km_t2_xa0']
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'

# # grid info
# g = xr.open_dataset(Ldir['grid'] / 'grid.nc')
# h = g.h.values
# h[g.mask_rho.values==0] = np.nan
# xrho = g.lon_rho.values
# yrho = g.lat_rho.values
# xp, yp = pfun.get_plon_plat(xrho,yrho)
# xu = g.lon_u.values
# yu = g.lat_u.values
# xv = g.lon_v.values
# yv = g.lat_v.values

fig, [ax1,ax2] = plt.subplots(1,2,figsize=(16,8))

# take a subset of the data so that we are averaging over an integer number of spring neap cycles at the end
# use 7 spring neap cycles, starting at index 257 and going to 2741 - these are the peaks in the 5km Qprism but similar for the other models
start_avg_ind = 257
end_avg_ind = 2741
start_avg_dt = '2020-09-13 05:00:00'
end_avg_dt = '2020-12-25 16:00:00' #this is tef_df.index[2740], the time indexing doesn't seem to be open half interval

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

    # Next, get the average salinity of the estuary and inner basin
    # These use the segment extractions for collections cei
    # the inner basin is segment b5i_p
    # for the whole estuary, need to get the weighted average of segment b5i_p and b5i_m
    seg_est_inner_ds = xr.open_dataset(seg_est_inner_fn)
    sbar_inner_ts = seg_est_inner_ds.salt.sel(seg='b5i_p',time=slice(start_avg_dt,end_avg_dt)).values
    sbar_outer_sill_ts = seg_est_inner_ds.salt.sel(seg='b5i_m',time=slice(start_avg_dt,end_avg_dt)).values
    vol_inner_ts = seg_est_inner_ds.volume.sel(seg='b5i_p',time=slice(start_avg_dt,end_avg_dt)).values
    vol_outer_sill_ts = seg_est_inner_ds.volume.sel(seg='b5i_m',time=slice(start_avg_dt,end_avg_dt)).values
    sbar_est_ts = ((sbar_inner_ts*vol_inner_ts)+(sbar_outer_sill_ts*vol_outer_sill_ts))/(vol_inner_ts+vol_outer_sill_ts)
    # get the time average values
    sbar_inner = sbar_inner_ts.mean()
    sbar_est = sbar_est_ts.mean()
    # to get the freshwater and saltwater fraction, define socn
    socn = 34
    #Calculate the freshwater and saltwater flushing times
    flushing_fresh_est_days[i] = ((socn-sbar_est)/(socn-sout_est))*flushing_est_days[i]
    flushing_fresh_inner_days[i] = ((socn-sbar_inner)/(socn-sout_inner))*flushing_inner_days[i]
    flushing_salt_est_days[i] = (sbar_est/sout_est)*flushing_est_days[i]
    flushing_salt_inner_days[i] = (sbar_inner)/(sout_inner)*flushing_inner_days[i]

print('Estuary flushing [days]: ',flushing_est_days)
print('Inner basin flushing [days]: ',flushing_inner_days)
print('Estuary freshwater flushing [days]: ',flushing_fresh_est_days)
print('Inner basin freshwater flushing [days]: ',flushing_fresh_inner_days)
print('Estuary saltwater flushing [days]: ',flushing_salt_est_days)
print('Inner basin saltwater flushing [days]: ',flushing_salt_inner_days)

#plot
ax1.plot(silllens_plot,flushing_est_days,c='k',marker='o',label='Volume flushing')
ax1.plot(silllens_plot,flushing_fresh_est_days,c='tab:cyan',marker='o',label='Freshwater flushing')
ax1.plot(silllens_plot,flushing_salt_est_days,c='tab:olive',marker='o',label='Saltwater flushing')
ax2.plot(silllens_plot,flushing_inner_days,c='k',marker='o',label='Volume flushing')
ax2.plot(silllens_plot,flushing_fresh_inner_days,c='tab:cyan',marker='o',label='Freshwater flushing')
ax2.plot(silllens_plot,flushing_salt_inner_days,c='tab:olive',marker='o',label='Saltwater flushing')

#add plot elements
ax1.set_xlabel('Sill length')
ax1.set_ylabel('Flushing time [days]')
ax1.set_ylim(0,85)
ax1.set_ylim(30,75)
ax1.set_title('Estuary')
ax1.grid(True)

ax2.set_xlabel('Sill length')
# ax2.set_ylabel('Flushing time [days]')
ax2.set_ylim(0,85)
ax2.set_ylim(30,75)
ax2.set_title('Inner basin')
ax2.grid(True)
ax2.legend()
plt.suptitle('Volume, freshwater, and saltwater flushing times')

# h, l = ax.get_legend_handles_labels()
# ph = [plt.plot([],marker="", ls="")[0]]*2
# handles = ph + h
# labels = [r'Inner basin reflux:', r'Outer basin reflux:'] + l
# # ax.legend(handles, labels, ncol=6)
# ax.legend(handles,labels,ncol=5)

fn_fig = Ldir['LOo'] / 'plots' / 'flushing_times.png' #UNCOMMENT TO PLOT
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

