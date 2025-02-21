"""
Compare particle based and TEF based efflux-reflux calculations
This version (2) uses an updated adjustment procedure following Hager et al 2022
This ensures that the salinities satisfy the inequalities in Cokelet and Stewart 1985
Also, here we only use the tidally averaged particles which is more representative of efflux-reflux


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
Q1_adj=np.zeros(5)
Q2_adj=np.zeros(5)
Q3_adj=np.zeros(5)
Q4_adj=np.zeros(5)
S1_adj=np.zeros(5)
S2_adj=np.zeros(5)
S3_adj=np.zeros(5)
S4_adj=np.zeros(5)
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

#first, satisfy salinity inequalities
#start by filling the adjusted variable arrays with the original TEF averages
Q1_adj[:]=Q1
Q2_adj[:]=Q2
Q3_adj[:]=Q3
Q4_adj[:]=Q4
S1_adj[:]=S1
S2_adj[:]=S2
S3_adj[:]=S3
S4_adj[:]=S4
#create arrays to store maximum and minimum allowable values for the salinities
#fill these when testing the four salinity inequalities
#we might be able to skip this if the conditions are met after the initial adjustment
# S1_minlim = np.nan*np.ones(4) #minimum limits on the salinities
# S2_minlim = np.nan*np.ones(4)
# S3_minlim = np.nan*np.ones(4)
# S2_maxlim = np.nan*np.ones(4) #maximum limits on the salinities
# S3_maxlim = np.nan*np.ones(4)
# S4_maxlim = np.nan*np.ones(4)

#adjust the salinities and transports to satisfy the salinity inequalities and volume and salt conservation
#since these are small arrays, we can use a for loop :)
for i in range(len(gctags)):

    #check the salinity inequalities, and if it is not met, set the salinities to the mean value of the two
    #first inequality: S4<=S2
    S2S4_mean = (S2_adj[i]+S4_adj[i])/2
    if S4_adj[i]>S2_adj[i]:
        print('S4>S2 for ',silllens[i],', adjusting!')
        S2_adj[i] = S2S4_mean + 0.0001
        S4_adj[i] = S2S4_mean - 0.0001
        Q2_adj[i] = (Q2[i]*S2[i])/S2_adj[i] #adjust volumes so that salt transports stay the same
        Q4_adj[i] = (Q4[i]*S4[i])/S4_adj[i]
        print('original values: S2 =',S2[i],'S4 =',S4[i],'Q2 =',Q2[i],'Q4 =',Q4[i])
        print('adjusted values: S2a=',S2_adj[i],'S4a=',S4_adj[i],'Q2a=',Q2_adj[i],'Q4a=',Q4_adj[i])

    #second inequality: S2<S1
    S1S2_mean = (S1_adj[i]+S2_adj[i])/2
    if S2_adj[i]>S1_adj[i]:
        print('S2>S1 for ',silllens[i],', adjusting!')
        S1_adj[i] = S1S2_mean + 0.0001
        S2_adj[i] = S1S2_mean - 0.0001
        Q1_adj[i] = (Q1[i]*S1[i])/S1_adj[i] #adjust volumes so that salt transports stay the same
        Q2_adj[i] = (Q2[i]*S2[i])/S2_adj[i]
        print('original values: S1 =',S1[i],'S2 =',S2[i],'Q1 =',Q1[i],'Q2 =',Q2[i])
        print('adjusted values: S1a=',S1_adj[i],'S2a=',S2_adj[i],'Q1a=',Q1_adj[i],'Q2a=',Q2_adj[i])  

    #third inequality: S4<S3
    S3S4_mean = (S3_adj[i]+S4_adj[i])/2
    if S4_adj[i]>S3_adj[i]:
        print('S4>S3 for ',silllens[i],', adjusting!')
        S3_adj[i] = S3S4_mean + 0.0001
        S4_adj[i] = S3S4_mean - 0.0001
        Q3_adj[i] = (Q3[i]*S3[i])/S3_adj[i] #adjust volumes so that salt transports stay the same
        Q4_adj[i] = (Q4[i]*S4[i])/S4_adj[i]
        print('original values: S3 =',S3[i],'S4 =',S4[i],'Q3 =',Q3[i],'Q4 =',Q4[i])
        print('adjusted values: S3a=',S3_adj[i],'S4a=',S4_adj[i],'Q3a=',Q3_adj[i],'Q4a=',Q4_adj[i])

    #fourth inequality: S3<=S1
    S1S3_mean = (S1_adj[i]+S3_adj[i])/2
    if S3_adj[i]>S1_adj[i]:
        print('S3>S1 for ',silllens[i],', adjusting!')
        S1_adj[i] = S1S3_mean + 0.0001
        S3_adj[i] = S1S3_mean - 0.0001
        Q1_adj[i] = (Q1[i]*S1[i])/S1_adj[i] #adjust volumes so that salt transports stay the same
        Q3_adj[i] = (Q3[i]*S3[i])/S3_adj[i]
        print('original values: S1 =',S1[i],'S3 =',S3[i],'Q1 =',Q1[i],'Q3 =',Q3[i])
        print('adjusted values: S1a=',S1_adj[i],'S3a=',S3_adj[i],'Q1a=',Q1_adj[i],'Q3a=',Q3_adj[i])

#next, check the volume conservation with the adjusted transports
vol_residual = Q1_adj-Q2_adj-Q3_adj+Q4_adj
#adjust volumes for perfect conservation
Q1_adj=Q1_adj-0.25*vol_residual #distribute the error evenly between the four transports, different signs are due to direction
Q2_adj=Q2_adj+0.25*vol_residual
Q3_adj=Q3_adj+0.25*vol_residual
Q4_adj=Q4_adj-0.25*vol_residual
print('After adjustment for perfect volume conservation:')
print('original values: Q1 =',Q1)
print('adjusted values: Q1a=',Q1_adj)
print('original values: Q2 =',Q2)
print('adjusted values: Q2a=',Q2_adj)
print('original values: Q3 =',Q3)
print('adjusted values: Q3a=',Q3_adj)
print('original values: Q4 =',Q4)
print('adjusted values: Q4a=',Q4_adj)

#next, check the salt conservation with the updated transports and salinities
salt_residual_adj = (Q1_adj*S1_adj)-(Q2_adj*S2_adj)-(Q3_adj*S3_adj)+(Q4_adj*S4_adj)
#allocate the salt transport error evenly between the four transports, and find the salinities necessary to keep the volume transports as Q_adj
S1_adj = (Q1_adj*S1_adj - 0.25*salt_residual_adj)/Q1_adj
S2_adj = (Q2_adj*S2_adj + 0.25*salt_residual_adj)/Q2_adj
S3_adj = (Q3_adj*S3_adj + 0.25*salt_residual_adj)/Q3_adj
S4_adj = (Q4_adj*S4_adj - 0.25*salt_residual_adj)/Q4_adj
print('After adjustment for perfect salt conservation:')
print('original values: S1 =',S1)
print('adjusted values: S1a=',S1_adj)
print('original values: S2 =',S2)
print('adjusted values: S2a=',S2_adj)
print('original values: S3 =',S3)
print('adjusted values: S3a=',S3_adj)
print('original values: S4 =',S4)
print('adjusted values: S4a=',S4_adj)

#check that the salinity inequalities are still satisfied
#if not, will need to add additional salt adjustment options to the code
print('S4<=S2 : ',S4<S2)
print('S2<S1 : ',S2<S1)
print('S4<S3 : ',S4<S3)
print('S3<=S1 : ',S3<S1)

# #check volume and salt conservation
# vol_residual = Q1-Q2-Q3+Q4
# salt_residual = (Q1*S1)-(Q2*S2)-(Q3*S3)+(Q4*S4)

# print('\nvolume residuals:')
# print(vol_residual)
# print('\nsalt residuals:')
# print(salt_residual)

# #adjust volumes for perfect conservation
# Q1_adj=Q1-0.25*vol_residual #distribute the error evenly between the four transports, different signs are due to direction
# Q2_adj=Q2+0.25*vol_residual
# Q3_adj=Q3+0.25*vol_residual
# Q4_adj=Q4-0.25*vol_residual

# #find the salt residual with the new adjusted volume transports
# salt_residual_Qadj = (Q1_adj*S1)-(Q2_adj*S2)-(Q3_adj*S3)+(Q4_adj*S4)

# #allocate the salt transport error evenly between the four transports, and find the salinities necessary to keep the volume transports as Q_adj
# S1_adj = (Q1_adj*S1 - 0.25*salt_residual_Qadj)/Q1_adj
# S2_adj = (Q2_adj*S2 + 0.25*salt_residual_Qadj)/Q2_adj
# S3_adj = (Q3_adj*S3 + 0.25*salt_residual_Qadj)/Q3_adj
# S4_adj = (Q4_adj*S4 - 0.25*salt_residual_Qadj)/Q4_adj

#calculate alphas with adjusted values
alpha_21_adj = (Q2_adj/Q1_adj)*((S2_adj-S4_adj)/(S1_adj-S4_adj))
alpha_31_adj = (Q3_adj/Q1_adj)*((S3_adj-S4_adj)/(S1_adj-S4_adj))
alpha_24_adj = (Q2_adj/Q4_adj)*((S1_adj-S2_adj)/(S1_adj-S4_adj))
alpha_34_adj = (Q3_adj/Q4_adj)*((S1_adj-S3_adj)/(S1_adj-S4_adj))

print('got alphas from tef')
print('After all adjustments:')
print('original values: alpha_21 =',alpha_21_basic)
print('adjusted values: alpha_21a=',alpha_21_adj)
print('original values: alpha_31 =',alpha_31_basic)
print('adjusted values: alpha_31a=',alpha_31_adj)
print('original values: alpha_24 =',alpha_24_basic)
print('adjusted values: alpha_24a=',alpha_24_adj)
print('original values: alpha_34 =',alpha_34_basic)
print('adjusted values: alpha_34a=',alpha_34_adj)

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

#Only plot the reflux fractions
ax.plot(silllens_plot,alpha_34_basic,marker='o',c='tab:pink',ls=':',label=r'TEF $\alpha_{34}$')
ax.plot(silllens_plot,alpha_21_basic,marker='o',c='tab:cyan',ls=':',label=r'TEF $\alpha_{21}$')
ax.plot(silllens_plot,alpha_34_adj,marker='o',c='tab:pink',ls='-',label=r'Adjusted TEF $\alpha_{34}$')
ax.plot(silllens_plot,alpha_21_adj,marker='o',c='tab:cyan',ls='-',label=r'Adjusted TEF $\alpha_{21}$')
ax.plot(silllens_plot,alpha_34_par_ta,marker='o',c='tab:pink',ls='--',label=r'Tidally averaged particles $\alpha_{34}$')
ax.plot(silllens_plot,alpha_21_par_ta,marker='o',c='tab:cyan',ls='--',label=r'Tidally averaged particles $\alpha_{21}$')

ax.set_xlabel('Sill length [km]')
ax.set_ylabel('Reflux fractions')
ax.set_xlim(0,85)
ax.set_ylim(-0.3,1.2)
ax.set_title('Efflux/reflux coefficients from TEF and particle tracking')
ax.grid(True)

h, l = ax.get_legend_handles_labels()
ph = [plt.plot([],marker="", ls="")[0]]*4
handles = ph + h
labels = [r'Inner basin reflux $\alpha_{34}$:', r'Outer basin reflux $\alpha_{21}$:',r'Efflux from inner basin $\alpha_{24}$:',r'Efflux from outer basin $\alpha_{31}$:'] + l
ax.legend(handles, labels, ncol=4, loc='lower right')

fn_fig = Ldir['LOo'] / 'plots' / 'efflux_reflux_coefs_compare_methods2.png' #UNCOMMENT TO PLOT
plt.savefig(fn_fig)
plt.close()

