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
gctags=['sill5km_c0', 'sill10km_c0', 'sill20kmdeep_c0', 'sill40km_c0', 'sill80km_c0']
gtagexs=['sill5km_t0_xa0', 'sill10km_t2_xa0', 'sill20kmdeep_t2_xa0', 'sill40km_t2_xa0', 'sill80km_t2_xa0']
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'

# grid info
g = xr.open_dataset(Ldir['grid'] / 'grid.nc')
h = g.h.values
h[g.mask_rho.values==0] = np.nan
xrho = g.lon_rho.values
yrho = g.lat_rho.values
xp, yp = pfun.get_plon_plat(xrho,yrho)
xu = g.lon_u.values
yu = g.lat_u.values
xv = g.lon_v.values
yv = g.lat_v.values

#Loop over sill lengths
for i in range(len(gctags)):
    #model and extraction info
    print(silllens[i])
    gctag=gctags[i]
    gtagex=gtagexs[i]
    in_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('bulk_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
    
    #Section b1
    sect_name = sect_1
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)         
    # adjust units
    tef_df['Q_p'] = tef_df['q_p']/1000
    tef_df['Q_m'] = tef_df['q_m']/1000
    tef_df['Q_prism']=tef_df['qprism']/1000
    # get variables for efflux-reflux calculation
    Q1 = tef_df['Q_p']
    Q2 = -tef_df['Q_m']
    S1 = tef_df['salt_p']
    S2 = tef_df['salt_m']

    #Section b5
    sect_name = sect_2
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
    # adjust units
    tef_df['Q_p'] = tef_df['q_p']/1000
    tef_df['Q_m'] = tef_df['q_m']/1000
    tef_df['Q_prism']=tef_df['qprism']/1000
    # get variables for efflux-reflux calculation
    Q3 = tef_df['Q_p']
    Q4 = -tef_df['Q_m']
    S3 = tef_df['salt_p']
    S4 = tef_df['salt_m']

    #Section b3 - we are using this to get an estimate of salinities in the upper and lower layers on the sill
    sect_name = sect_mid
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
    # adjust units
    tef_df['Q_p'] = tef_df['q_p']/1000
    tef_df['Q_m'] = tef_df['q_m']/1000
    tef_df['Q_prism']=tef_df['qprism']/1000
    # get variables for storage term in efflux-reflux calculation
    S_bottom = tef_df['salt_p']
    S_top = tef_df['salt_m']

    #Get the volume of the sill area
    #get grid stuff
    #find the volume between the two longitudes
    #divide by 2

    # #check volume and salt conservation
    # vol_residual = Q1-Q2-Q3+Q4
    # salt_residual = (Q1*S1)-(Q2*S2)-(Q3*S3)+(Q4*S4)

    #calculate alphas basic without adjustment or storage term
    alpha_21_basic = (Q2/Q1)*((S2-S4)/(S1-S4))
    alpha_31_basic = (Q3/Q1)*((S3-S4)/(S1-S4))
    alpha_24_basic = (Q2/Q4)*((S1-S2)/(S1-S4))
    alpha_34_basic = (Q3/Q4)*((S1-S3)/(S1-S4))

    #calculate storage term with centered differences #make sure to use 3600s and volume divided by 1000



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

