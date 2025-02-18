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

Qin_A_mean=np.zeros(5)
Qout_A_mean=np.zeros(5)
sin_A_mean=np.zeros(5)
sout_A_mean=np.zeros(5)
Qin_B_mean=np.zeros(5)
Qout_B_mean=np.zeros(5)
sin_B_mean=np.zeros(5)
sout_B_mean=np.zeros(5)
Q1=np.zeros(5)
Q2=np.zeros(5)
Q3=np.zeros(5)
Q4=np.zeros(5)
S1=np.zeros(5)
S2=np.zeros(5)
S3=np.zeros(5)
S4=np.zeros(5)
alpha_24_basic=np.zeros(5)
alpha_34_basic=np.zeros(5)
alpha_31_basic=np.zeros(5)
alpha_21_basic=np.zeros(5)
silllens_plot=[5,10,20,40,80]
alpha_34_td_top_avg=np.zeros(5)
alpha_21_td_top_avg=np.zeros(5)
alpha_34_td_bottom_avg=np.zeros(5)
alpha_21_td_bottom_avg=np.zeros(5)
alpha_34_basic_avg=np.zeros(5)
alpha_21_basic_avg=np.zeros(5)

sect_1 ='b1'
sect_2='b5'
sect_mid='b3'

plot_color = ['tab:red','tab:orange','tab:green','tab:blue','tab:purple'] #COLORS FOR 5 models

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

fig, ax = plt.subplots(1,1,figsize=(15,8))

# take a subset of the data so that we are averaging over an integer number of spring neap cycles at the end
# use 7 spring neap cycles, starting at index 257 and going to 2741 - these are the peaks in the 5km Qprism but similar for the other models
start_avg_ind = 257
end_avg_ind = 2741

#parameters for smoothing
sg_window_size = 13
sg_order = 5

#Loop over sill lengths
for i in range(len(gctags)):
# for i in range(len(gctags)-1):
    #model and extraction info
    print(silllens[i])
    gctag=gctags[i]
    gtagex=gtagexs[i]
    gridname=gridnames[i]
    in_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('bulk_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
    
    #Section b1
    sect_name = sect_1
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)         
    # adjust units
    # tef_df['Q_p'] = tef_df['q_p']/1000
    # tef_df['Q_m'] = tef_df['q_m']/1000
    # tef_df['Q_prism']=tef_df['qprism']/1000
    # get variables for efflux-reflux calculation
    Q1 = tef_df['q_p'][start_avg_ind:end_avg_ind] #keep original units for comparing with V_sill #maybe should add .values???
    Q2 = -tef_df['q_m'][start_avg_ind:end_avg_ind]
    S1 = tef_df['salt_p'][start_avg_ind:end_avg_ind]
    S2 = tef_df['salt_m'][start_avg_ind:end_avg_ind]

    Q1_smooth = scipy.signal.savgol_filter(tef_df['q_p'],sg_window_size,sg_order)[start_avg_ind:end_avg_ind]
    Q2_smooth = scipy.signal.savgol_filter(-tef_df['q_m'],sg_window_size,sg_order)[start_avg_ind:end_avg_ind]
    S1_smooth = scipy.signal.savgol_filter(tef_df['salt_p'],sg_window_size,sg_order)[start_avg_ind:end_avg_ind]
    S2_smooth = scipy.signal.savgol_filter(tef_df['salt_m'],sg_window_size,sg_order)[start_avg_ind:end_avg_ind]

    #Section b5
    sect_name = sect_2
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
    # adjust units
    # tef_df['Q_p'] = tef_df['q_p']/1000
    # tef_df['Q_m'] = tef_df['q_m']/1000
    # tef_df['Q_prism']=tef_df['qprism']/1000
    # get variables for efflux-reflux calculation
    Q3 = tef_df['q_p'][start_avg_ind:end_avg_ind] #keep original units for comparing with V_sill
    Q4 = -tef_df['q_m'][start_avg_ind:end_avg_ind]
    S3 = tef_df['salt_p'][start_avg_ind:end_avg_ind]
    S4 = tef_df['salt_m'][start_avg_ind:end_avg_ind]

    Q3_smooth = scipy.signal.savgol_filter(tef_df['q_p'],sg_window_size,sg_order)[start_avg_ind:end_avg_ind]
    Q4_smooth = scipy.signal.savgol_filter(-tef_df['q_m'],sg_window_size,sg_order)[start_avg_ind:end_avg_ind]
    S3_smooth = scipy.signal.savgol_filter(tef_df['salt_p'],sg_window_size,sg_order)[start_avg_ind:end_avg_ind]
    S4_smooth = scipy.signal.savgol_filter(tef_df['salt_m'],sg_window_size,sg_order)[start_avg_ind:end_avg_ind]

    #Section b3 - we are using this to get an estimate of salinities in the upper and lower layers on the sill
    sect_name = sect_mid
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
    # # adjust units
    # tef_df['Q_p'] = tef_df['q_p']/1000
    # tef_df['Q_m'] = tef_df['q_m']/1000
    # tef_df['Q_prism']=tef_df['qprism']/1000
    # get variables for storage term in efflux-reflux calculation
    S_bottom = tef_df['salt_p'][start_avg_ind-1:end_avg_ind+1] #keep one extra value on each side for the salt
    S_top = tef_df['salt_m'][start_avg_ind-1:end_avg_ind+1] #this way when we do centered differences the length will be the same as the other variables

    #first, clip out where salinity changes more than 0.5 psu in an hour using a median filter
    #if a value is more than 0.5psu different from the median with the two surrounding points, replace with the average of the two nearest points
    #NEW: if a value is more than 0.5psu different from 5 point median, replace with the median
    S_bottom_values = tef_df['salt_p'].values
    S_top_values = tef_df['salt_m'].values
    # S_bottom_clip = np.where(np.abs(scipy.signal.medfilt(S_bottom_values)-S_bottom_values)>0.5,np.concatenate(([np.nan],(S_bottom_values[:-2]+S_bottom_values[2:])/2,[np.nan])),S_bottom_values)
    # S_top_clip = np.where(np.abs(scipy.signal.medfilt(S_top_values)-S_top_values)>0.5,np.concatenate(([np.nan],(S_top_values[:-2]+S_top_values[2:])/2,[np.nan])),S_top_values)
    S_bottom_clip = np.where(np.abs(scipy.signal.medfilt(S_bottom_values,5)-S_bottom_values)>0.5,np.signal.medfilt(S_bottom_values,5),S_bottom_values)
    S_top_clip = np.where(np.abs(scipy.signal.medfilt(S_top_values,5)-S_top_values)>0.5,np.signal.medfilt(S_top_values,5),S_top_values)
    #now, smooth with savgol filter
    # S_bottom_smooth = scipy.signal.savgol_filter(tef_df['salt_p'],sg_window_size,sg_order)[start_avg_ind-1:end_avg_ind+1]
    # S_top_smooth = scipy.signal.savgol_filter(tef_df['salt_m'],sg_window_size,sg_order)[start_avg_ind-1:end_avg_ind+1]
    S_bottom_smooth = scipy.signal.savgol_filter(S_bottom_clip,sg_window_size,sg_order)[start_avg_ind-1:end_avg_ind+1]
    S_top_smooth = scipy.signal.savgol_filter(S_top_clip,sg_window_size,sg_order)[start_avg_ind-1:end_avg_ind+1]

    # dSdt_bottom_smooth = (1/3600)*scipy.signal.savgol_filter(tef_df['salt_p'],sg_window_size,sg_order,deriv=1)[start_avg_ind:end_avg_ind] #use the filter to take the derivative
    # dSdt_top_smooth = (1/3600)*scipy.signal.savgol_filter(tef_df['salt_m'],sg_window_size,sg_order,deriv=1)[start_avg_ind:end_avg_ind]
    dSdt_bottom_smooth = (1/3600)*scipy.signal.savgol_filter(S_bottom_clip,sg_window_size,sg_order,deriv=1)[start_avg_ind:end_avg_ind] #use the filter to take the derivative
    dSdt_top_smooth = (1/3600)*scipy.signal.savgol_filter(S_top_clip,sg_window_size,sg_order,deriv=1)[start_avg_ind:end_avg_ind]

    #Get the volume of the sill area
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
    #get the b1 part of the sect_df
    sect_df_b1=sect_df[sect_df.sn=='b1']
    #get the index of the u values
    b1_ind_u=sect_df_b1.i.values[0]
    #get the longitude of the section from lon_u
    b1_lon_u=lon_u[0,b1_ind_u]
    #do the same for b5
    sect_df_b5=sect_df[sect_df.sn=='b5']
    b5_ind_u=sect_df_b5.i.values[0]
    b5_lon_u=lon_u[0,b5_ind_u]
    #mask the h array ouside of the sill area
    h_sill = g.h.values
    h_sill[mask_rho==0] = np.nan #mask the land
    h_sill[lon_rho<b1_lon_u] = np.nan #mask seaward of sill
    h_sill[lon_rho>b5_lon_u] = np.nan #mask landward of sill
    #get the sill volume
    V_sill = np.nansum(DA*h_sill)
    #divide by 2
    V_top = V_sill/2 #for now, assume layer volume is half of sill
    V_bottom = V_sill/2
    print(V_sill)

    # #check volume and salt conservation
    # vol_residual = Q1-Q2-Q3+Q4
    # salt_residual = (Q1*S1)-(Q2*S2)-(Q3*S3)+(Q4*S4)

    #calculate alphas basic from timeseries (not means) without adjustment or storage term
    alpha_21_basic_timeseries = (Q2/Q1)*((S2-S4)/(S1-S4))
    alpha_31_basic_timeseries = (Q3/Q1)*((S3-S4)/(S1-S4))
    alpha_24_basic_timeseries = (Q2/Q4)*((S1-S2)/(S1-S4))
    alpha_34_basic_timeseries = (Q3/Q4)*((S1-S3)/(S1-S4))

    #calculate storage term with centered differences #make sure to use 3600s and volume divided by 1000
    ddt_S_top = (S_top.values[2:]-S_top.values[:-2])/(2*3600)
    # ddt_S_top = np.concatenate(([np.nan],ddt_S_top,[np.nan]))
    storage_21_top = (1/(Q1*(S1-S4)))*V_top*ddt_S_top
    storage_24_top = (1/(Q4*(S4-S1)))*V_top*ddt_S_top

    #calculate time dependent alphas
    alpha_21_td_top = alpha_21_basic_timeseries + storage_21_top
    alpha_24_td_top = alpha_24_basic_timeseries + storage_24_top
    alpha_31_td_top = 1-alpha_21_td_top
    alpha_34_td_top = 1-alpha_24_td_top

    #to test the method, try calculating the same terms from the bottom layer budget
    ddt_S_bottom = (S_bottom.values[2:]-S_bottom.values[:-2])/(2*3600)
    storage_31_bottom = (1/(Q1*(S1-S4)))*V_bottom*ddt_S_bottom    
    storage_34_bottom = (1/(Q4*(S4-S1)))*V_bottom*ddt_S_bottom

    alpha_31_td_bottom = alpha_31_basic_timeseries + storage_31_bottom
    alpha_34_td_bottom = alpha_34_basic_timeseries + storage_34_bottom
    alpha_21_td_bottom = 1-alpha_31_td_bottom
    alpha_24_td_bottom = 1-alpha_34_td_bottom

    #get averages for summary plot
    alpha_34_td_top_avg[i] = alpha_34_td_top.mean()
    alpha_21_td_top_avg[i] = alpha_21_td_top.mean()
    alpha_34_td_bottom_avg[i] = alpha_34_td_bottom.mean()
    alpha_21_td_bottom_avg[i] = alpha_21_td_bottom.mean()

    #get basic (unadjusted) time average value
    alpha_21_basic_avg[i] = (Q2.mean()/Q1.mean())*((S2.mean()-S4.mean())/(S1.mean()-S4.mean()))
    # alpha_31_basic_avg[i] = (Q3.mean()/Q1.mean())*((S3.mean()-S4.mean())/(S1.mean()-S4.mean()))
    # alpha_24_basic_avg[i] = (Q2.mean()/Q4.mean())*((S1.mean()-S2.mean())/(S1.mean()-S4.mean()))
    alpha_34_basic_avg[i] = (Q3.mean()/Q4.mean())*((S1.mean()-S3.mean())/(S1.mean()-S4.mean()))

    #count times when salinity condition is not met
    print('S4>S3 count: ',(S4>S3).sum()) #technically some of these should be inequalities
    print('S3>S1 count: ',(S3>S1).sum())
    print('S4>S2 count: ',(S4>S2).sum())
    print('S2>S1 count: ',(S2>S1).sum())

    #print some values
    print('outer reflux w/o storage')
    print(alpha_21_basic_timeseries.mean())
    print('outer reflux w/ storage (top layer calc)')
    print(alpha_21_td_top.mean())
    print('outer reflux w/ storage (bottom layer calc)')
    print(alpha_21_td_bottom.mean())
    print('inner reflux w/o storage')
    print(alpha_34_basic_timeseries.mean())
    print('inner reflux w/ storage (top layer calc)')
    print(alpha_34_td_top.mean())
    print('inner reflux w/ storage (bottom layer calc)')
    print(alpha_34_td_bottom.mean())

    #plot
    plot_time = tef_df.index[start_avg_ind:end_avg_ind]
    # ax.plot(plot_time,alpha_34_td_top,ls='-',c=plot_color[i],label=r'$\alpha_{34}$ '+silllens[i])
    # ax.plot(plot_time,alpha_21_td_top,ls='--',c=plot_color[i],label=r'$\alpha_{21}$ '+silllens[i])
    ax.plot(plot_time,alpha_34_td_bottom,ls='-',c=plot_color[i],label=r'$\alpha_{34}$ '+silllens[i])
    ax.plot(plot_time,alpha_21_td_bottom,ls='--',c=plot_color[i],label=r'$\alpha_{21}$ '+silllens[i])

    #try getting smoothed version
    #calculate alphas basic from timeseries (not means) without adjustment or storage term
    alpha_21_basic_timeseries_smooth = (Q2_smooth/Q1_smooth)*((S2_smooth-S4_smooth)/(S1_smooth-S4_smooth))
    alpha_31_basic_timeseries_smooth = (Q3_smooth/Q1_smooth)*((S3_smooth-S4_smooth)/(S1_smooth-S4_smooth))
    alpha_24_basic_timeseries_smooth = (Q2_smooth/Q4_smooth)*((S1_smooth-S2_smooth)/(S1_smooth-S4_smooth))
    alpha_34_basic_timeseries_smooth = (Q3_smooth/Q4_smooth)*((S1_smooth-S3_smooth)/(S1_smooth-S4_smooth))

    #calculate storage term with centered differences #make sure to use 3600s
    ddt_S_top_smooth_cd = (S_top_smooth[2:]-S_top_smooth[:-2])/(2*3600) #get derivative from centered differences on smoothed data, can also use filter dervative
    # ddt_S_top = np.concatenate(([np.nan],ddt_S_top,[np.nan]))
    storage_21_top_smooth = (1/(Q1_smooth*(S1_smooth-S4_smooth)))*V_top*dSdt_top_smooth
    storage_24_top_smooth = (1/(Q4_smooth*(S4_smooth-S1_smooth)))*V_top*dSdt_top_smooth

    #calculate time dependent alphas
    alpha_21_td_top_smooth = alpha_21_basic_timeseries_smooth + storage_21_top_smooth
    alpha_24_td_top_smooth = alpha_24_basic_timeseries_smooth + storage_24_top_smooth
    alpha_31_td_top_smooth = 1-alpha_21_td_top_smooth
    alpha_34_td_top_smooth = 1-alpha_24_td_top_smooth

    #to test the method, try calculating the same terms from the bottom layer budget
    ddt_S_bottom_smooth_cd = (S_bottom_smooth[2:]-S_bottom_smooth[:-2])/(2*3600)
    storage_31_bottom_smooth = (1/(Q1_smooth*(S1_smooth-S4_smooth)))*V_bottom*dSdt_bottom_smooth    
    storage_34_bottom_smooth = (1/(Q4_smooth*(S4_smooth-S1_smooth)))*V_bottom*dSdt_bottom_smooth

    alpha_31_td_bottom_smooth = alpha_31_basic_timeseries_smooth + storage_31_bottom_smooth
    alpha_34_td_bottom_smooth = alpha_34_basic_timeseries_smooth + storage_34_bottom_smooth
    alpha_21_td_bottom_smooth = 1-alpha_31_td_bottom_smooth
    alpha_24_td_bottom_smooth = 1-alpha_34_td_bottom_smooth


#add plot elements
ax.set_xlabel('Time')
ax.set_ylabel('Reflux coefficient')
# ax.set_ylim(-4,4)
ax.set_ylim(0,1)
# ax.set_title('Time-dependent efflux/reflux coefficients')
ax.set_title('Time-dependent efflux/reflux coefficients from bottom layer budget')
ax.grid(True)
ax.legend(ncol=5)
ax.set_xlim('2020-09-01','2021-01-01')

h, l = ax.get_legend_handles_labels()
ph = [plt.plot([],marker="", ls="")[0]]*2
handles = ph + h
labels = [r'Inner basin reflux:', r'Outer basin reflux:'] + l
# ax.legend(handles, labels, ncol=6)
ax.legend(handles,labels,ncol=5)

fn_fig = Ldir['LOo'] / 'plots' / 'efflux_reflux_coefs_time_dependent.png' #UNCOMMENT TO PLOT
plt.savefig(fn_fig)
plt.close()

#The problem is in the storage_21 term, but what part?
# fig, [ax1,ax2,ax3,ax4,ax5,ax6] = plt.subplots(6,1,figsize=(16,16))
# ax1.plot(plot_time,ddt_S_top,ls='-',label=r'$d/dt(S_{top})$')
# ax2.plot(plot_time,S_top[1:-1],ls='-',label=r'$S_{top}$')
# ax3.plot(plot_time,Q1,ls='-',label=r'Q1')
# # ax3.plot(plot_time,1/Q1,ls='-',label=r'1/Q1')
# ax4.plot(plot_time,S1-S4,ls='-',label=r'S1-S4')
# ax5.plot(plot_time,S1,label=r'S1')
# ax6.plot(plot_time,S4,label=r'S4')
# # ax5.plot(plot_time,1/(S1-S4),ls='-',label=r'1/(S1-S4)')
# # ax6.plot(plot_time,storage_21_top,ls='-',label=r'storage_21_top')
# ax1.legend()
# ax2.legend()
# ax3.legend()
# ax4.legend()
# ax5.legend()
# ax6.legend()
# ax1.set_ylim(-0.0001,0.0001)
# # ax6.set_ylim(-5,5)
# plt.suptitle('Components of top layer alpha_21 storage term')

fig, axs = plt.subplots(9,2,figsize=(16,16),sharey='row',sharex=True)
axs[0,0].plot(plot_time,alpha_21_td_top,ls='-',label=r'Time dependent \alpha_{21}')
axs[1,0].plot(plot_time,alpha_21_basic_timeseries,ls='-',label=r'Basic \alpha_{21} timeseries')
axs[2,0].plot(plot_time,storage_21_top,ls='-',label=r'storage term')
axs[3,0].plot(plot_time,ddt_S_top,ls='-',label=r'$d/dt(S_{top})$')
axs[4,0].plot(plot_time,S_top[1:-1],ls='-',label=r'$S_{top}$')
axs[5,0].plot(plot_time,Q1,ls='-',label=r'Q1')
# ax6.plot(plot_time,1/Q1,ls='-',label=r'1/Q1')
axs[6,0].plot(plot_time,S1-S4,ls='-',label=r'S1-S4')
axs[7,0].plot(plot_time,S1,label=r'S1')
axs[8,0].plot(plot_time,S4,label=r'S4')

axs[0,1].plot(plot_time,alpha_21_td_top_smooth,ls='-',label=r'Time dependent \alpha_{21}')
axs[1,1].plot(plot_time,alpha_21_basic_timeseries_smooth,ls='-',label=r'Basic \alpha_{21} timeseries')
axs[2,1].plot(plot_time,storage_21_top_smooth,ls='-',label=r'storage term')
axs[3,1].plot(plot_time,ddt_S_top_smooth_cd,ls='-',label=r'$d/dt(S_{top})$ smoothed + centered differences')
axs[3,1].plot(plot_time,dSdt_top_smooth,ls='-',label=r'$d/dt(S_{top})$ derivative from SG filter')
axs[4,1].plot(plot_time,S_top_smooth[1:-1],ls='-',label=r'$S_{top}$')
axs[5,1].plot(plot_time,Q1_smooth,ls='-',label=r'Q1')
# ax3.plot(plot_time,1/Q1,ls='-',label=r'1/Q1')
axs[6,1].plot(plot_time,S1_smooth-S4_smooth,ls='-',label=r'S1-S4')
axs[7,1].plot(plot_time,S1_smooth,label=r'S1')
axs[8,1].plot(plot_time,S4_smooth,label=r'S4')
# ax5.plot(plot_time,1/(S1-S4),ls='-',label=r'1/(S1-S4)')
# ax6.plot(plot_time,storage_21_top,ls='-',label=r'storage_21_top')
axs[0,0].legend()
axs[1,0].legend()
axs[2,0].legend()
axs[3,0].legend()
axs[4,0].legend()
axs[5,0].legend()
axs[6,0].legend()
axs[7,0].legend()
axs[8,0].legend()
axs[0,1].legend()
axs[1,1].legend()
axs[2,1].legend()
axs[3,1].legend()
axs[4,1].legend()
axs[5,1].legend()
axs[6,1].legend()
axs[7,1].legend()
axs[8,1].legend()
# axs[0,0].set_ylim(-0.0001,0.0001)
# ax6.set_ylim(-5,5)
plt.suptitle(r'Components of top layer time dependent \alpha_{21}')
axs[0,0].set_title('Raw data')
axs[0,1].set_title('Savitzky-Golay filter (window length:'+str(sg_window_size)+', order:'+str(sg_order))

fn_fig = Ldir['LOo'] / 'plots' / 'efflux_reflux_coefs_time_dependent_storagetest.png' #UNCOMMENT TO PLOT
plt.savefig(fn_fig)
plt.close()

#Averages plot
fig, [ax1,ax2] = plt.subplots(1,2,figsize=(16,8))
#plot
ax1.plot(silllens_plot,alpha_21_td_top_avg,c='tab:orange',marker='o',label='Time-dependent (from top layer budget)')
ax1.plot(silllens_plot,alpha_21_td_bottom_avg,c='tab:blue',marker='o',label='Time-dependent (from bottom layer budget)')
ax1.plot(silllens_plot,alpha_21_basic_avg,c='tab:grey',marker='o',ls='--',label='Time-averaged')
ax2.plot(silllens_plot,alpha_34_td_top_avg,c='tab:orange',marker='o')
ax2.plot(silllens_plot,alpha_34_td_bottom_avg,c='tab:blue',marker='o')
ax2.plot(silllens_plot,alpha_34_basic_avg,c='tab:grey',marker='o',ls='--')

#add plot elements
ax1.set_xlabel('Sill length')
ax1.set_ylabel('Average reflux coefficient')
ax1.set_xlim(0,85)
ax1.set_ylim(-0.1,0.9)
ax1.set_title(r'Outer basin reflux $\alpha_{21}$')
ax1.grid(True)
ax1.legend()

ax2.set_xlabel('Sill length')
# ax2.set_ylabel('Flushing time [days]')
ax2.set_xlim(0,85)
ax2.set_ylim(-0.1,0.9)
ax2.set_title(r'Inner basin reflux $\alpha_{34}$')
ax2.grid(True)
# ax2.legend()
plt.suptitle('Time-dependent reflux coefficients')

fn_fig = Ldir['LOo'] / 'plots' / 'efflux_reflux_coefs_time_dependent_avg.png' #UNCOMMENT TO PLOT
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

