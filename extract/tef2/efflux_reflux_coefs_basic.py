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

# sect_list = [item.name for item in in_dir.glob('*.p')]
# if Ldir['testing']:
#     sect_list = ['ss2.p']
#sect_list = ['a3.p','b3.p','c3.p']
#sect_list = ['a1.p','a2.p','a3.p','a4.p','a5.p','b1.p','b2.p','b3.p','b4.p','b5.p','c1.p','c2.p','c3.p','c4.p','c5.p']
#sect_list = ['a1.p','a3.p','b1.p','b2.p','b3.p','b4.p','b5.p','c3.p']
#plot_label = ['a1','a3','b1','b2','b3','b4','b5','c3']
#plot_color = ['tab:red','tab:orange','tab:olive','tab:green','tab:cyan','tab:blue','tab:purple','tab:pink']
#plot_color = ['k','tab:gray','tab:red','tab:orange','tab:green','tab:cyan','tab:blue','tab:brown']

#sect_list = ['b1','b2','b3','b4','b5']
# sect_choice = 'b5'
# sect_choice = 'b3'
sect_1 ='b1'
sect_2='b5'

#plot_label = ['b1','b2','b3','b4','b5']
#plot_color = ['tab:red','tab:orange','tab:green','tab:cyan','tab:blue']
plot_color = ['tab:red','tab:orange','tab:green','tab:blue','tab:purple'] #COLORS FOR 5 models
# plot_color = ['tab:red','tab:green','tab:purple']

# silllens=['5km','20km','80km']
# gctags=['sill5km_c0','sill20kmdeep_c0','sill80km_c0']
# gtagexs=['sill5km_t0_xa0', 'sill20kmdeep_t2_xa0', 'sill80km_t2_xa0']
silllens=['5km', '10km', '20km', '40km', '80km']
gctags=['sill5km_c0', 'sill10km_c0', 'sill20kmdeep_c0', 'sill40km_c0', 'sill80km_c0']
gtagexs=['sill5km_t0_xa0', 'sill10km_t2_xa0', 'sill20kmdeep_t2_xa0', 'sill40km_t2_xa0', 'sill80km_t2_xa0']
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'

#in_dir = out_dir0 / ('bulk_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
# ds01s = ['2020.09.01_2020.12.31','2020.09.15_2020.11.15','2020.09.15_2020.11.15']
# out_dir = out_dir0 / ('bulk_plots_multimodel_' + Ldir['ds0'] + '_' + Ldir['ds1'])
# Lfun.make_dir(out_dir, clean=True)

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

# GET AVERAGE TEF VARIABLES FROM BOTH SECTIONS

#Section A (b1)
for i in range(len(gctags)):
    print(silllens[i])
    gctag=gctags[i]
    gtagex=gtagexs[i]
    # in_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('bulk_hourly_' + ds01s[i])
    in_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('bulk_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
    #sect_name = sect_list[i]
    sect_name = sect_1
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))

    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
            
    # adjust units
    tef_df['Q_p'] = tef_df['q_p']/1000
    tef_df['Q_m'] = tef_df['q_m']/1000
    tef_df['Q_prism']=tef_df['qprism']/1000

    peak_list, peak_props = scipy.signal.find_peaks(tef_df['Q_prism'])
    print('First and last qprism peaks:')
    print(peak_list[1])
    print(peak_list[-1])
    print(peak_list[-1]-peak_list[1])

    # get tide info from the tide excursion calculator
    excur_dir = out_dir0 / ('tide_excursion_' + Ldir['ds0'] + '_' + Ldir['ds1'])
    te_fn = excur_dir / ('TE_b3.p') #could change if using other sections
    TE = pd.read_pickle(te_fn)

    # Qindeltas=tef_df['Q_p']*(tef_df['salt_p']-tef_df['salt_m'])
    # Qindeltas_avg = Qindeltas.mean()
    # Qindeltas_spring = Qindeltas.loc[TE['t_spring']]
    # Qindeltas_neap = Qindeltas.loc[TE['t_neap']]
    # print('\n')
    # print(gctag)
    # print('\nAverage Qindeltas:')
    # print(Qindeltas_avg)
    # print('\nSpring Qindeltas:')
    # print(Qindeltas_spring)
    # print('\nNeap Qindeltas:')
    # print(Qindeltas_neap)

    Qin_A_mean[i] = tef_df['Q_p'].mean() #CHANGE THESE TO USE INTEGER SN CYCLES
    Qout_A_mean[i] = tef_df['Q_m'].mean()
    sin_A_mean[i] = tef_df['salt_p'].mean()
    sout_A_mean[i] = tef_df['salt_m'].mean()

    Q1[i] = np.abs(tef_df['Q_p'].mean()) #ADD ABSOLUTE VALUE TO WORK WITH THE FORMULAS (NEED TO CHECK THAT THIS DOESN'T LOSE SIGN INFO)
    Q2[i] = np.abs(tef_df['Q_m'].mean())
    S1[i] = tef_df['salt_p'].mean()
    S2[i] = tef_df['salt_m'].mean()

    # print('\nQ_in_A (Q1):')
    # print(Qin_A_mean[i])
    # print('\nQ_out_A (Q2):')
    # print(Qout_A_mean[i])
    # print('\ns_in_A (S1):')
    # print(sin_A_mean[i])
    # print('\ns_out_A (S2):')
    # print(sout_A_mean[i])
                    
#Section 2 (b5)
for i in range(len(gctags)):
    # print(silllens[i])
    gctag=gctags[i]
    gtagex=gtagexs[i]
    # in_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('bulk_hourly_' + ds01s[i])
    in_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('bulk_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
    #sect_name = sect_list[i]
    sect_name = sect_2
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))

    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
            
    # adjust units
    tef_df['Q_p'] = tef_df['q_p']/1000
    tef_df['Q_m'] = tef_df['q_m']/1000
    tef_df['Q_prism']=tef_df['qprism']/1000

    peak_list, peak_props = scipy.signal.find_peaks(tef_df['Q_prism'])
    print('First and last qprism peaks:')
    print(peak_list[1])
    print(peak_list[-1])
    print(peak_list[-1]-peak_list[1])

    # get tide info from the tide excursion calculator
    excur_dir = out_dir0 / ('tide_excursion_' + Ldir['ds0'] + '_' + Ldir['ds1'])
    te_fn = excur_dir / ('TE_b3.p') #could change if using other sections
    TE = pd.read_pickle(te_fn)

    # Qindeltas=tef_df['Q_p']*(tef_df['salt_p']-tef_df['salt_m'])
    # Qindeltas_avg = Qindeltas.mean()
    # Qindeltas_spring = Qindeltas.loc[TE['t_spring']]
    # Qindeltas_neap = Qindeltas.loc[TE['t_neap']]
    # print('\n')
    # print(gctag)
    # print('\nAverage Qindeltas:')
    # print(Qindeltas_avg)
    # print('\nSpring Qindeltas:')
    # print(Qindeltas_spring)
    # print('\nNeap Qindeltas:')
    # print(Qindeltas_neap)

    Qin_B_mean[i] = tef_df['Q_p'].mean()
    Qout_B_mean[i] = tef_df['Q_m'].mean()
    sin_B_mean[i] = tef_df['salt_p'].mean()
    sout_B_mean[i] = tef_df['salt_m'].mean()

    Q3[i] = np.abs(tef_df['Q_p'].mean())
    Q4[i] = np.abs(tef_df['Q_m'].mean())
    S3[i] = tef_df['salt_p'].mean()
    S4[i] = tef_df['salt_m'].mean()

    # print('\nQ_in_B (Q3):')
    # print(Qin_B_mean[i])
    # print('\nQ_out_B (Q4):')
    # print(Qout_B_mean[i])
    # print('\ns_in_B (S3):')
    # print(sin_B_mean[i])
    # print('\ns_out_B (S4):')
    # print(sout_B_mean[i])

#check volume and salt conservation
vol_residual = Q1-Q2-Q3+Q4
salt_residual = (Q1*S1)-(Q2*S2)-(Q3*S3)+(Q4*S4)

# print('\nvolume residuals:')
# print(vol_residual)
# print('\nsalt residuals:')
# print(salt_residual)

#calculate alphas
alpha_21_basic = (Q2/Q1)*((S2-S4)/(S1-S4))
alpha_31_basic = (Q3/Q1)*((S3-S4)/(S1-S4))
alpha_24_basic = (Q2/Q4)*((S1-S2)/(S1-S4))
alpha_34_basic = (Q3/Q4)*((S1-S3)/(S1-S4))

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

