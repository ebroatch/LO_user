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

#Loop over sill lengths
for i in range(len(gctags)):
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
    Q1 = tef_df['q_p'] #keep original units for comparing with V_sill
    Q2 = -tef_df['q_m']
    S1 = tef_df['salt_p']
    S2 = tef_df['salt_m']

    #Section b5
    sect_name = sect_2
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
    # adjust units
    # tef_df['Q_p'] = tef_df['q_p']/1000
    # tef_df['Q_m'] = tef_df['q_m']/1000
    # tef_df['Q_prism']=tef_df['qprism']/1000
    # get variables for efflux-reflux calculation
    Q3 = tef_df['q_p'] #keep original units for comparing with V_sill
    Q4 = -tef_df['q_m']
    S3 = tef_df['salt_p']
    S4 = tef_df['salt_m']

    #Section b3 - we are using this to get an estimate of salinities in the upper and lower layers on the sill
    sect_name = sect_mid
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
    # # adjust units
    # tef_df['Q_p'] = tef_df['q_p']/1000
    # tef_df['Q_m'] = tef_df['q_m']/1000
    # tef_df['Q_prism']=tef_df['qprism']/1000
    # get variables for storage term in efflux-reflux calculation
    S_bottom = tef_df['salt_p']
    S_top = tef_df['salt_m']

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

    #calculate alphas basic without adjustment or storage term
    alpha_21_basic = (Q2/Q1)*((S2-S4)/(S1-S4))
    alpha_31_basic = (Q3/Q1)*((S3-S4)/(S1-S4))
    alpha_24_basic = (Q2/Q4)*((S1-S2)/(S1-S4))
    alpha_34_basic = (Q3/Q4)*((S1-S3)/(S1-S4))

    #calculate storage term with centered differences #make sure to use 3600s and volume divided by 1000
    ddt_S_top = (S_top.values[2:]-S_top.values[:-2])/(2*3600)
    ddt_S_top = np.concatenate(([np.nan],ddt_S_top,[np.nan]))
    storage_21 = (1/(Q1*(S1-S4)))*V_top*ddt_S_top
    storage_24 = (1/(Q4*(S4-S1)))*V_top*ddt_S_top

    #calculate time dependent alphas
    alpha_21_td = alpha_21_basic + storage_21
    alpha_24_td = alpha_24_basic + storage_24
    alpha_34_td = 1-alpha_24_td

    #plot
    ax.plot(tef_df.index,alpha_34_td,ls='-',c=plot_color[i],label=r'Inner basin reflux $\alpha_{34} '+silllens[i])
    ax.plot(tef_df.index,alpha_21_td,ls='--',c=plot_color[i],label=r'Outer basin reflux $\alpha_{21} '+silllens[i])

#add plot elements
ax.set_xlabel('Time')
ax.set_ylabel('Reflux coefficient')
ax.set_ylim(-10,10)
ax.set_title('Time-dependent efflux/reflux coefficients')
ax.grid(True)
ax.legend()

fn_fig = Ldir['LOo'] / 'plots' / 'efflux_reflux_coefs_time_dependent.png' #UNCOMMENT TO PLOT
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

