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
# use 7 spring neap cycles, starting at index 257 and going to 2741 - these are the peaks in the 5km Qprism but similar for the other models
start_avg_ind = 257
end_avg_ind = 2741

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
    # print('First and last qprism peaks:')
    # print(peak_list[1])
    # print(peak_list[-1])
    # print(peak_list[-1]-peak_list[1])

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

    Qin_A_mean[i] = tef_df['Q_p'].mean() #This is the basic mean of whole timeseries
    Qout_A_mean[i] = tef_df['Q_m'].mean()
    sin_A_mean[i] = tef_df['salt_p'].mean()
    sout_A_mean[i] = tef_df['salt_m'].mean()

    Q1[i] = tef_df['Q_p'][start_avg_ind:end_avg_ind].mean() #average over 7 s/n cycles
    Q2[i] = -tef_df['Q_m'][start_avg_ind:end_avg_ind].mean() #average over 7 s/n cycles and change sign so that all Q are positive
    S1[i] = tef_df['salt_p'][start_avg_ind:end_avg_ind].mean()
    S2[i] = tef_df['salt_m'][start_avg_ind:end_avg_ind].mean()

    # print('\nQ_in_A, Q1):')
    # print(Qin_A_mean[i], Q1[i])
    # print('\nQ_out_A, Q2:')
    # print(Qout_A_mean[i], Q2[i])
    # print('\ns_in_A, S1:')
    # print(sin_A_mean[i], S1[i])
    # print('\ns_out_A, S2:')
    # print(sout_A_mean[i], S2[i])
                    
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
    # print('First and last qprism peaks:')
    # print(peak_list[1])
    # print(peak_list[-1])
    # print(peak_list[-1]-peak_list[1])

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

    Q3[i] = tef_df['Q_p'][start_avg_ind:end_avg_ind].mean() #average over 7 s/n cycles
    Q4[i] = -tef_df['Q_m'][start_avg_ind:end_avg_ind].mean() #average over 7 s/n cycles and change sign so that all Q are positive
    S3[i] = tef_df['salt_p'][start_avg_ind:end_avg_ind].mean()
    S4[i] = tef_df['salt_m'][start_avg_ind:end_avg_ind].mean()

    # print('\nQ_in_B, Q3:')
    # print(Qin_B_mean[i], Q3[i])
    # print('\nQ_out_B, Q4:')
    # print(Qout_B_mean[i], Q4[i])
    # print('\ns_in_B, S3:')
    # print(sin_B_mean[i], S3[i])
    # print('\ns_out_B, S4:')
    # print(sout_B_mean[i], S4[i])

#calculate basic alphas with no adjustment
alpha_21_basic = (Q2/Q1)*((S2-S4)/(S1-S4))
alpha_31_basic = (Q3/Q1)*((S3-S4)/(S1-S4))
alpha_24_basic = (Q2/Q4)*((S1-S2)/(S1-S4))
alpha_34_basic = (Q3/Q4)*((S1-S3)/(S1-S4))

print('\nalpha_21 (outer basin reflux): ')
print(alpha_21_basic)
print('\nalpha_31 (efflux from outer basin): ')
print(alpha_31_basic)
print('\nalpha_34 (inner basin reflux): ')
print(alpha_34_basic)
print('\nalpha_24 (efflux from inner basin): ')
print(alpha_24_basic)

#check volume and salt conservation
vol_residual = Q1-Q2-Q3+Q4
salt_residual = (Q1*S1)-(Q2*S2)-(Q3*S3)+(Q4*S4)

# print('\nvolume residuals:')
# print(vol_residual)
# print('\nsalt residuals:')
# print(salt_residual)

#adjust volumes for perfect conservation
Q1_adj=Q1-0.25*vol_residual #distribute the error evenly between the four transports, different signs are due to direction
Q2_adj=Q2+0.25*vol_residual
Q3_adj=Q3+0.25*vol_residual
Q4_adj=Q4-0.25*vol_residual

#find the salt residual with the new adjusted volume transports
salt_residual_Qadj = (Q1_adj*S1)-(Q2_adj*S2)-(Q3_adj*S3)+(Q4_adj*S4)

#allocate the salt transport error evenly between the four transports, and find the salinities necessary to keep the volume transports as Q_adj
S1_adj = (Q1_adj*S1 - 0.25*salt_residual_Qadj)/Q1_adj
S2_adj = (Q2_adj*S2 + 0.25*salt_residual_Qadj)/Q2_adj
S3_adj = (Q3_adj*S3 + 0.25*salt_residual_Qadj)/Q3_adj
S4_adj = (Q4_adj*S4 - 0.25*salt_residual_Qadj)/Q4_adj

#calculate alphas with adjusted values
alpha_21_adj = (Q2_adj/Q1_adj)*((S2_adj-S4_adj)/(S1_adj-S4_adj))
alpha_31_adj = (Q3_adj/Q1_adj)*((S3_adj-S4_adj)/(S1_adj-S4_adj))
alpha_24_adj = (Q2_adj/Q4_adj)*((S1_adj-S2_adj)/(S1_adj-S4_adj))
alpha_34_adj = (Q3_adj/Q4_adj)*((S1_adj-S3_adj)/(S1_adj-S4_adj))

print('\nalpha_21 (adjusted outer basin reflux): ')
print(alpha_21_adj)
print('\nalpha_31 (adjusted efflux from outer basin): ')
print(alpha_31_adj)
print('\nalpha_34 (adjusted inner basin reflux): ')
print(alpha_34_adj)
print('\nalpha_24 (adjusted efflux from inner basin): ')
print(alpha_24_adj)

print('\nalpha_21+alpha_31 (adjusted outer basin fractions): ')
print(alpha_21_adj+alpha_31_adj)
print('\nalpha_24+alpha_34 (adjusted outer basin fractions): ')
print(alpha_24_adj+alpha_34_adj)

print('\nvolume discrepancy due to adjustment [%]: ')
print('Q1: ',100*(Q1_adj-Q1)/Q1)
print('Q2: ',100*(Q2_adj-Q2)/Q2)
print('Q3: ',100*(Q3_adj-Q3)/Q3)
print('Q4: ',100*(Q4_adj-Q4)/Q4)

print('\nsalinity discrepancy due to adjustment [%]: ')
print('S1: ',100*(S1_adj-S1)/S1)
print('S2: ',100*(S2_adj-S2)/S2)
print('S3: ',100*(S3_adj-S3)/S3)
print('S4: ',100*(S4_adj-S4)/S4)

print('\nsalt transport discrepancy due to adjustment [%]: ')
print('Q1*S1: ',100*(Q1_adj*S1_adj-Q1*S1)/Q1*S1)
print('Q2*S2: ',100*(Q2_adj*S2_adj-Q2*S2)/Q2*S2)
print('Q3*S3: ',100*(Q3_adj*S3_adj-Q3*S3)/Q3*S3)
print('Q4*S4: ',100*(Q4_adj*S4_adj-Q4*S4)/Q4*S4)

print('\nalpha discrepancy due to adjustment [%]: ')
print('alpha_21: ',100*(alpha_21_adj-alpha_21_basic)/alpha_21_basic)
print('alpha_31: ',100*(alpha_31_adj-alpha_31_basic)/alpha_31_basic)
print('alpha_24: ',100*(alpha_24_adj-alpha_24_basic)/alpha_24_basic)
print('alpha_34: ',100*(alpha_34_adj-alpha_34_basic)/alpha_34_basic)

fig, ax = plt.subplots(1,1,figsize=(15,8))

ax.plot(silllens_plot,alpha_34_adj,marker='o',c='tab:pink',ls='-',label=r'Inner basin reflux $\alpha_{34}$')
ax.plot(silllens_plot,alpha_21_adj,marker='o',c='tab:cyan',ls='-',label=r'Outer basin reflux $\alpha_{21}$')
ax.plot(silllens_plot,alpha_24_adj,marker='o',c=plt.cm.tab20(13),ls='--',label=r'Efflux from inner basin $\alpha_{24}$')
ax.plot(silllens_plot,alpha_31_adj,marker='o',c=plt.cm.tab20(19),ls='--',label=r'Efflux from outer basin $\alpha_{31}$')

ax.set_xlabel('Sill length [km]')
ax.set_ylabel('Efflux/reflux coefficients')
ax.set_xlim(0,85)
ax.set_ylim(-0.1,1.1)
ax.set_title('Efflux/reflux fractions from TEF')
ax.grid(True)
ax.legend()

fn_fig = Ldir['LOo'] / 'plots' / 'efflux_reflux_coefs_TEF_adjust.png' #UNCOMMENT TO PLOT
plt.savefig(fn_fig)
plt.close()

