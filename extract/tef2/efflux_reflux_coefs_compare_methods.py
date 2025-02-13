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
from lo_user_tools import llxyfun
import tef_fun
import scipy.signal

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']


########## EFFLUX-REFLUX FROM TEF ##########
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

silllens=['5km', '10km', '20km', '40km', '80km']
gctags=['sill5km_c0', 'sill10km_c0', 'sill20kmdeep_c0', 'sill40km_c0', 'sill80km_c0']
gtagexs=['sill5km_t0_xa0', 'sill10km_t2_xa0', 'sill20kmdeep_t2_xa0', 'sill40km_t2_xa0', 'sill80km_t2_xa0']
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'

# GET AVERAGE TEF VARIABLES FROM BOTH SECTIONS
# use 7 spring neap cycles, starting at index 257 and going to 2741 - these are the peaks in the 5km Qprism but similar for the other models
start_avg_ind = 257
end_avg_ind = 2741

#Section A (b1)
for i in range(len(gctags)):
    gctag=gctags[i]
    gtagex=gtagexs[i]
    in_dir = Ldir['LOo'] / 'extract' / gtagex / 'tef2' / ('bulk_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
    sect_name = sect_1
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))

    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
            
    # adjust units
    tef_df['Q_p'] = tef_df['q_p']/1000
    tef_df['Q_m'] = tef_df['q_m']/1000
    tef_df['Q_prism']=tef_df['qprism']/1000

    peak_list, peak_props = scipy.signal.find_peaks(tef_df['Q_prism'])

    Q1[i] = tef_df['Q_p'][start_avg_ind:end_avg_ind].mean() #average over 7 s/n cycles
    Q2[i] = -tef_df['Q_m'][start_avg_ind:end_avg_ind].mean() #average over 7 s/n cycles and change sign so that all Q are positive
    S1[i] = tef_df['salt_p'][start_avg_ind:end_avg_ind].mean()
    S2[i] = tef_df['salt_m'][start_avg_ind:end_avg_ind].mean()

                    
#Section 2 (b5)
for i in range(len(gctags)):
    gctag=gctags[i]
    gtagex=gtagexs[i]
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

    Q3[i] = tef_df['Q_p'][start_avg_ind:end_avg_ind].mean() #average over 7 s/n cycles
    Q4[i] = -tef_df['Q_m'][start_avg_ind:end_avg_ind].mean() #average over 7 s/n cycles and change sign so that all Q are positive
    S3[i] = tef_df['salt_p'][start_avg_ind:end_avg_ind].mean()
    S4[i] = tef_df['salt_m'][start_avg_ind:end_avg_ind].mean()


#calculate basic alphas with no adjustment
alpha_21_basic = (Q2/Q1)*((S2-S4)/(S1-S4))
alpha_31_basic = (Q3/Q1)*((S3-S4)/(S1-S4))
alpha_24_basic = (Q2/Q4)*((S1-S2)/(S1-S4))
alpha_34_basic = (Q3/Q4)*((S1-S3)/(S1-S4))

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

print('got alphas from tef')

########## PARTICLES (NOT TIDALLY AVERAGED) ##########
alpha_24_par=np.zeros(5)
alpha_34_par=np.zeros(5)
alpha_31_par=np.zeros(5)
alpha_21_par=np.zeros(5)

for i in range(5):
    #loop over runs
    if i==0:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(45e3,0,45)
        linecolor = 'tab:red'
        silllenlabel = '5km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill5km_t0_xa0/sill5kmest_3d/release_2020.09.01.nc'
        tef_5km_fn = Ldir['LOo'] / 'extract/sill5km_t0_xa0/tef2/bulk_hourly_2020.09.01_2020.12.31' #this is for the Qprism timekeeper and shading
    elif i==1:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(50e3,0,45)
        linecolor = 'tab:orange'
        silllenlabel = '10km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill10km_t2_xa0/sill10kmest_3d/release_2020.09.01.nc'
    elif i==2:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(60e3,0,45)
        linecolor = 'tab:green'
        silllenlabel = '20km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill20kmdeep_t2_xa0/sill20kmdeepest_3d/release_2020.09.01.nc'
    elif i==3:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(80e3,0,45)
        linecolor = 'tab:blue'
        silllenlabel = '40km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill40km_t2_xa0/sill40kmest_3d/release_2020.09.01.nc'
    elif i==4:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(120e3,0,45)
        linecolor = 'tab:purple'
        silllenlabel = '80km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill80km_t2_xa0/sill80kmest_3d/release_2020.09.01.nc'

    # get Datasets
    dsr = xr.open_dataset(fn, decode_times=False)
    NT, NP = dsr.lon.shape
    #longitude data
    lon_vals = dsr.lon.values #longitudes of the particles
    #time data
    time_hours = dsr.Time.values
    dsr.close()

    #find when particles are in different regions
    lon_in = lon_vals >= sillland #boolean array of particles in the inner basin over time
    lon_sill = (lon_vals > sillsea) & (lon_vals < sillland) #boolean array of particles on the sill over time

    #make array with code number for each region
    region_codes = (3*lon_in.astype(int))+(2*lon_sill.astype(int)) #NEW CODES: this gives 0 for outer basin and ocean, 2 for sill, and 3 for inner basin
    #find transitions between regions
    region_codes_transition = np.diff(region_codes,axis=0) #WITH NEW CODES this gives 0 staying in same region,-1 inner to sill,+1 sill to inner,+2 outer to sill,-2 sill to outer,+3 outer to inner direct,-3 inner to outer direct
    #get the transitions all in a row with zeros removed
    region_codes_transition_nan = np.where(region_codes_transition==0,np.nan,region_codes_transition) #change 0 to nan so that we can remove them and only look at consecutive transitions
    a = (~np.isnan(region_codes_transition_nan)).argsort(0, kind='mergesort') #this gives the indices to sort the array with the nans first along the time axis, use mergesort to preserve order of other elements
    region_codes_transition_consecutive = region_codes_transition_nan[a, np.arange(a.shape[1])[None,:]] #this should sort all the nans to the top of the column

    #now we need to search for different patterns within the columns which indicate visits to the sill as efflux or reflux, and count them
    #NEW CODES:
    #-1,-2 is in->sill->out (efflux)
    #-1,+1 is in->sill->in (reflux to inner basin)
    #+2,+1 is out->sill->in (efflux)
    #+2,-2 is out->sill->out (reflux to outer basin)
    #-3 is direct in->out (transit sill in <1h)
    #+3 is direct out->in (transit sill in <1h)
    pattern_inout = [-1,-2]
    pattern_inin = [-1,1]
    pattern_outin = [2,1]
    pattern_outout = [2,-2]
    direct_inout = -3
    direct_outin = +3

    inout_bool = (region_codes_transition_consecutive[:-1,:]==pattern_inout[0]) & (region_codes_transition_consecutive[1:,:]==pattern_inout[1]) #boolean arrays of where the patterns are found
    inin_bool = (region_codes_transition_consecutive[:-1,:]==pattern_inin[0]) & (region_codes_transition_consecutive[1:,:]==pattern_inin[1])
    outin_bool = (region_codes_transition_consecutive[:-1,:]==pattern_outin[0]) & (region_codes_transition_consecutive[1:,:]==pattern_outin[1])
    outout_bool = (region_codes_transition_consecutive[:-1,:]==pattern_outout[0]) & (region_codes_transition_consecutive[1:,:]==pattern_outout[1])

    inout_count = np.sum(inout_bool)
    inin_count = np.sum(inin_bool)
    outin_count = np.sum(outin_bool)
    outout_count = np.sum(outout_bool)

    direct_inout_count = np.sum(region_codes_transition_consecutive==direct_inout) #these are included with the regular count when calculating the alphas (only matters for the 5km model)
    direct_outin_count = np.sum(region_codes_transition_consecutive==direct_outin)

    #next, find the efflux reflux fractions
    alpha_24 = (inout_count+direct_inout_count)/(inout_count+direct_inout_count+inin_count)
    alpha_34 = inin_count/(inout_count+direct_inout_count+inin_count)
    alpha_31 = (outin_count+direct_outin_count)/(outin_count+direct_outin_count+outout_count)
    alpha_21 = outout_count/(outin_count+direct_outin_count+outout_count) 

    alpha_24_par[i]=alpha_24
    alpha_34_par[i]=alpha_34
    alpha_31_par[i]=alpha_31
    alpha_21_par[i]=alpha_21

print('got alphas from particles')

########## PARTICLES (TIDALLY AVERAGED) ##########
alpha_24_par_ta=np.zeros(5)
alpha_34_par_ta=np.zeros(5)
alpha_31_par_ta=np.zeros(5)
alpha_21_par_ta=np.zeros(5)

for i in range(5):
    #loop over runs
    if i==0:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(45e3,0,45)
        linecolor = 'tab:red'
        silllenlabel = '5km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill5km_t0_xa0/sill5kmest_3d/release_2020.09.01.nc'
        tef_5km_fn = Ldir['LOo'] / 'extract/sill5km_t0_xa0/tef2/bulk_hourly_2020.09.01_2020.12.31' #this is for the Qprism timekeeper and shading
    elif i==1:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(50e3,0,45)
        linecolor = 'tab:orange'
        silllenlabel = '10km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill10km_t2_xa0/sill10kmest_3d/release_2020.09.01.nc'
    elif i==2:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(60e3,0,45)
        linecolor = 'tab:green'
        silllenlabel = '20km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill20kmdeep_t2_xa0/sill20kmdeepest_3d/release_2020.09.01.nc'
    elif i==3:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(80e3,0,45)
        linecolor = 'tab:blue'
        silllenlabel = '40km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill40km_t2_xa0/sill40kmest_3d/release_2020.09.01.nc'
    elif i==4:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(120e3,0,45)
        linecolor = 'tab:purple'
        silllenlabel = '80km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill80km_t2_xa0/sill80kmest_3d/release_2020.09.01.nc'

    # get Datasets
    dsr = xr.open_dataset(fn, decode_times=False)
    NT, NP = dsr.lon.shape
    #longitude data
    lon_vals = dsr.lon.values #longitudes of the particles
    lon_vals_ta = zfun.lowpass(lon_vals, f='godin')
    #time data
    time_hours = dsr.Time.values
    dsr.close()

    #find when particles are in different regions
    lon_in = lon_vals_ta >= sillland #boolean array of particles in the inner basin over time
    lon_sill = (lon_vals_ta > sillsea) & (lon_vals_ta < sillland) #boolean array of particles on the sill over time

    #make array with code number for each region
    region_codes = (3*lon_in.astype(int))+(2*lon_sill.astype(int)) #NEW CODES: this gives 0 for outer basin and ocean, 2 for sill, and 3 for inner basin
    #find transitions between regions
    region_codes_transition = np.diff(region_codes,axis=0) #WITH NEW CODES this gives 0 staying in same region,-1 inner to sill,+1 sill to inner,+2 outer to sill,-2 sill to outer,+3 outer to inner direct,-3 inner to outer direct
    #get the transitions all in a row with zeros removed
    region_codes_transition_nan = np.where(region_codes_transition==0,np.nan,region_codes_transition) #change 0 to nan so that we can remove them and only look at consecutive transitions
    a = (~np.isnan(region_codes_transition_nan)).argsort(0, kind='mergesort') #this gives the indices to sort the array with the nans first along the time axis, use mergesort to preserve order of other elements
    region_codes_transition_consecutive = region_codes_transition_nan[a, np.arange(a.shape[1])[None,:]] #this should sort all the nans to the top of the column

    #now we need to search for different patterns within the columns which indicate visits to the sill as efflux or reflux, and count them
    #NEW CODES:
    #-1,-2 is in->sill->out (efflux)
    #-1,+1 is in->sill->in (reflux to inner basin)
    #+2,+1 is out->sill->in (efflux)
    #+2,-2 is out->sill->out (reflux to outer basin)
    #-3 is direct in->out (transit sill in <1h)
    #+3 is direct out->in (transit sill in <1h)
    pattern_inout = [-1,-2]
    pattern_inin = [-1,1]
    pattern_outin = [2,1]
    pattern_outout = [2,-2]
    direct_inout = -3
    direct_outin = +3

    inout_bool = (region_codes_transition_consecutive[:-1,:]==pattern_inout[0]) & (region_codes_transition_consecutive[1:,:]==pattern_inout[1]) #boolean arrays of where the patterns are found
    inin_bool = (region_codes_transition_consecutive[:-1,:]==pattern_inin[0]) & (region_codes_transition_consecutive[1:,:]==pattern_inin[1])
    outin_bool = (region_codes_transition_consecutive[:-1,:]==pattern_outin[0]) & (region_codes_transition_consecutive[1:,:]==pattern_outin[1])
    outout_bool = (region_codes_transition_consecutive[:-1,:]==pattern_outout[0]) & (region_codes_transition_consecutive[1:,:]==pattern_outout[1])

    inout_count = np.sum(inout_bool)
    inin_count = np.sum(inin_bool)
    outin_count = np.sum(outin_bool)
    outout_count = np.sum(outout_bool)

    direct_inout_count = np.sum(region_codes_transition_consecutive==direct_inout) #these are included with the regular count when calculating the alphas (only matters for the 5km model)
    direct_outin_count = np.sum(region_codes_transition_consecutive==direct_outin)

    #next, find the efflux reflux fractions
    alpha_24 = (inout_count+direct_inout_count)/(inout_count+direct_inout_count+inin_count)
    alpha_34 = inin_count/(inout_count+direct_inout_count+inin_count)
    alpha_31 = (outin_count+direct_outin_count)/(outin_count+direct_outin_count+outout_count)
    alpha_21 = outout_count/(outin_count+direct_outin_count+outout_count) 

    alpha_24_par_ta[i]=alpha_24
    alpha_34_par_ta[i]=alpha_34
    alpha_31_par_ta[i]=alpha_31
    alpha_21_par_ta[i]=alpha_21

print('got alphas from tidally averaged particles')

########## PLOTTING ##########
fig, ax = plt.subplots(1,1,figsize=(15,8))

ax.plot(silllens_plot,alpha_34_adj,marker='o',c='tab:pink',ls='-',label=r'Adjusted TEF $\alpha_{34}$')
ax.plot(silllens_plot,alpha_21_adj,marker='o',c='tab:cyan',ls='-',label=r'Adjusted TEF $\alpha_{21}$')
ax.plot(silllens_plot,alpha_24_adj,marker='o',c=plt.cm.tab20(13),ls='-',label=r'Adjusted TEF $\alpha_{24}$')
ax.plot(silllens_plot,alpha_31_adj,marker='o',c=plt.cm.tab20(19),ls='-',label=r'Adjusted TEF $\alpha_{31}$')
ax.plot(silllens_plot,alpha_34_par_ta,marker='o',c='tab:pink',ls='--',label=r'Tidally averaged particles $\alpha_{34}$')
ax.plot(silllens_plot,alpha_21_par_ta,marker='o',c='tab:cyan',ls='--',label=r'Tidally averaged particles $\alpha_{21}$')
ax.plot(silllens_plot,alpha_24_par_ta,marker='o',c=plt.cm.tab20(13),ls='--',label=r'Tidally averaged particles $\alpha_{24}$')
ax.plot(silllens_plot,alpha_31_par_ta,marker='o',c=plt.cm.tab20(19),ls='--',label=r'Tidally averaged particles $\alpha_{31}$')
ax.plot(silllens_plot,alpha_34_par,marker='o',c='tab:pink',ls=':',label=r'Particles $\alpha_{34}$')
ax.plot(silllens_plot,alpha_21_par,marker='o',c='tab:cyan',ls=':',label=r'Particles $\alpha_{21}$')
ax.plot(silllens_plot,alpha_24_par,marker='o',c=plt.cm.tab20(13),ls=':',label=r'Particles $\alpha_{24}$')
ax.plot(silllens_plot,alpha_31_par,marker='o',c=plt.cm.tab20(19),ls=':',label=r'Particles $\alpha_{31}$')

ax.set_xlabel('Sill length [km]')
ax.set_ylabel('Efflux/reflux coefficients')
ax.set_xlim(0,85)
ax.set_ylim(-0.2,1.2)
ax.set_title('Efflux/reflux fractions')
ax.grid(True)

h, l = ax.get_legend_handles_labels()
ph = [plt.plot([],marker="", ls="")[0]]*4
handles = ph + h
labels = [r'Inner basin reflux $\alpha_{34}$:', r'Outer basin reflux $\alpha_{21}$:',r'Efflux from inner basin $\alpha_{24}$:',r'Efflux from outer basin $\alpha_{31}$:'] + l
ax.legend(handles, labels, ncol=4)

fn_fig = Ldir['LOo'] / 'plots' / 'efflux_reflux_coefs_compare_methods.png' #UNCOMMENT TO PLOT
plt.savefig(fn_fig)
plt.close()

