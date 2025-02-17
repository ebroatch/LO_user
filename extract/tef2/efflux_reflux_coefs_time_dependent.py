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

#Loop over sill lengths
# for i in range(len(gctags)):
for i in range(len(gctags)-1):
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

#Averages plot
fig, [ax1,ax2] = plt.subplots(1,2,figsize=(16,8))
#plot
ax1.plot(silllens_plot,alpha_21_td_top_avg,c='tab:orange',marker='o',label='Time-dependent (from top layer budget)')
ax1.plot(silllens_plot,alpha_21_td_bottom_avg,c='tab:blue',marker='o',label='Time-dependent (from bottom layer budget)')
ax1.plot(silllens_plot,alpha_21_basic_avg,c='tab:grey',marker='o',label='Time-averaged')
ax2.plot(silllens_plot,alpha_34_td_top_avg,c='tab:orange',marker='o')
ax2.plot(silllens_plot,alpha_34_td_bottom_avg,c='tab:blue',marker='o')
ax2.plot(silllens_plot,alpha_34_basic_avg,c='tab:grey',marker='o')

#add plot elements
ax1.set_xlabel('Sill length')
ax1.set_ylabel('Average reflux coefficient')
ax1.set_xlim(0,85)
ax1.set_ylim(0,1)
ax1.set_title(r'Outer basin reflux $\alpha_{21}$')
ax1.grid(True)
ax1.legend()

ax2.set_xlabel('Sill length')
# ax2.set_ylabel('Flushing time [days]')
ax2.set_xlim(0,85)
ax2.set_ylim(0,1)
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

