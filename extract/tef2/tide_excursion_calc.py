"""
Process tef2 extractions, giving transport vs. salinity for:
volume, salt, and other variables.

Can be run by making interactive choices in the terminal, or using command line
arguments (as when it is run by extract_sections.py).

PERFORMANCE: Using the new "binned_statistic" method this is significantly faster,
about 5 minutes for 83 section for a year with only salt. Is this really faster?

To test on mac:
run process_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.06 -test True
run process_sections -gtx cas6_v00Stock_uu0mb -ctag c0 -0 2021.07.04 -1 2021.07.06

And for a full year:

(this only has salt)
run process_sections -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31 -test True

"""

import sys
import xarray as xr
import numpy as np
import pickle
from time import time
import pandas as pd
from scipy.stats import binned_statistic
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from cmocean import cm

from lo_tools import plotting_functions as pfun

from lo_tools import Lfun, zrfun, zfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing
# use the arg -sect_name to pick section for plots

# import tef_fun

#choose closest date to save peaks from
neardayinput = input("Enter nearest date to use (YYYY-MM-DD): ")
nearday = np.datetime64(neardayinput)

#get extractions and processed sections, and set directories to save to
gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)
c_dir = Ldir['LOo'] / 'extract' / 'tef2' / ('sections' + '_' + gctag)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('tide_excursion_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

sect_list = [item.name for item in in_dir.glob('*.nc')]
if Ldir['testing']:
    #sect_list = ['b1.nc','b2.nc','b3.nc','b4.nc','b5.nc']
    sect_list=['b3.nc']
#sect_list = [Ldir['sect_name']+'.nc'] #section to plot #not sure if this syntax is correct
    
# make vn_list by inspecting the first section
ds = xr.open_dataset(in_dir / sect_list[0])
vn_list = [item for item in ds.data_vars \
    if (len(ds[item].dims) == 3) and (item not in ['vel','DZ'])]
ds.close()

print('\nCalculating parameters:')
#print(str(in_dir))

tt00 = time()
pad=36

#constants
beta=7.7e-4
g=9.81
socn=34
Cd=2.5e-3
omega=1.4e-4
QR=1000
T_M2=12.42 #period in hours
T_S2=12
T=(T_M2+T_S2)/2
Ts=T*3600 #period in seconds

#spring and neap times
# if Ldir['gridname']=='sill5km':
#     t_spring_ebb = pd.Timestamp('2020-10-12 16:00:00') #5km model
# else:
#     t_spring_ebb = pd.Timestamp('2020-10-12 17:00:00') #20km and 80km model

# if Ldir['gridname']=='sill80km':
#     t_spring_flood = pd.Timestamp('2020-10-12 12:00:00') #80km model
#     t_neap_flood = pd.Timestamp('2020-10-05 18:00:00') #80km model
# else:
#     t_spring_flood = pd.Timestamp('2020-10-12 10:00:00') #20km and 5km model
#     t_neap_flood = pd.Timestamp('2020-10-05 04:00:00') #20km and 5km model

# t_neap_ebb = pd.Timestamp('2020-10-05 11:00:00')

#list for tidal excursion
TE_spring=[]
TE_neap=[]
sect_tick=[]

#start parameter space plot
plt.close('all')
fs = 14
pfun.start_plot(fs=fs, figsize=(15,15))
fig, ax = plt.subplots(1, 1)

#loop over sections
for ext_fn in sect_list:
    tt0 = time()
    print(ext_fn)
    sect_label = ext_fn.replace('.nc','')

    # load extraction fields
    ds = xr.open_dataset(in_dir / ext_fn)

    # get variables from section extraction
    ot = ds['time'].to_numpy() #shape NT
    zeta = ds['zeta'].to_numpy() #shape NT,NX
    h = ds['h'].to_numpy() #shape NX
    dd = ds['dd'].to_numpy() #shape NX
    dz = ds['DZ'].to_numpy() #shape NT,NZ,NX
    u = ds['vel'].to_numpy() #shape NT,NZ,NX
    s = ds['salt'].to_numpy() #shape NT,NZ,NX

    # calculate derived variables
    dA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() #shape NT,NZ,NX
    A = np.sum(dA, axis=(1,2)) #shape NT
    q = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['vel'].to_numpy() #shape NT,NZ,NX
    qnet = np.sum(q, axis=(1,2)) #shape NT
    unet = qnet/A #shape NT

    # calculate section variables (not time varying)
    H_thalweg=np.max(h) #thalweg depth (use as depth for N0 and M)
    A_sect=np.sum(h*dd) #section area (not time varying)

    #calculate qprism
    qprism=zfun.lowpass(np.abs(qnet - zfun.lowpass(qnet, f='godin')), f='godin')/2

    #cut off nans for peak finding - scipy find_peaks can't handle nans
    qprism_pf=(qprism[pad:-pad+1])[pad:-pad+1]
    qnet_pf=(qnet[pad:-pad+1])[pad:-pad+1]
    unet_pf=(unet[pad:-pad+1])[pad:-pad+1]
    ot_pf=(ot[pad:-pad+1])[pad:-pad+1]

    #find peaks
    [spring_peaks, spring_properties]=find_peaks(qprism_pf) 
    [neap_peaks, neap_properties]=find_peaks(-qprism_pf)
    [flood_peaks, flood_properties]=find_peaks(qnet_pf, height=0) #need to use height kwarg (only include peaks above 0) so that the properties will include heights of peaks
    [ebb_peaks, ebb_properties]=find_peaks(-qnet_pf, height=0)

    #locate closest spring and neap to nearday, and closest flood and ebb to spring and neap
    spring_times = ot_pf[spring_peaks]
    neap_times = ot_pf[neap_peaks]
    t_spring = spring_times[np.argmin(np.abs(spring_times-nearday))]
    t_neap = neap_times[np.argmin(np.abs(neap_times-nearday))]
    
    flood_times = ot_pf[flood_peaks]
    ebb_times = ot_pf[ebb_peaks]
    t_spring_flood = flood_times[np.argmin(np.abs(flood_times-t_spring))]
    t_spring_ebb = ebb_times[np.argmin(np.abs(ebb_times-t_spring))]
    t_neap_flood = flood_times[np.argmin(np.abs(flood_times-t_neap))]
    t_neap_ebb = ebb_times[np.argmin(np.abs(ebb_times-t_neap))]

    #FINDING BIGGEST UT NEAR CHOSEN SPRING TIDE
    #change to use unet and not qnet
    [u_flood_peaks, u_flood_properties]=find_peaks(unet_pf, height=0) #need to use height kwarg (only include peaks above 0) so that the properties will include heights of peaks
    [u_ebb_peaks, u_ebb_properties]=find_peaks(-unet_pf, height=0)
    u_flood_times = ot_pf[u_flood_peaks]
    u_ebb_times = ot_pf[u_ebb_peaks]
    u_flood_heights=u_flood_properties['peak_heights']
    u_ebb_heights=u_ebb_properties['peak_heights']
    if u_flood_times[0]<u_ebb_times[0]:
        if len(u_flood_times)==len(u_ebb_times):
            UT1=(u_flood_heights+u_ebb_heights)/2 #UT1 is flood first pairs
            UT2=(u_flood_heights[1:]+u_ebb_heights[:-1])/2 #UT2 is ebb first pairs
            t_UT1=u_flood_times+((u_ebb_times-u_flood_times)/2) #average time of UT1 flood and ebb
            t_UT2=u_flood_times[1:]+((u_ebb_times[:-1]-u_flood_times[1:])/2) #average time of UT2 ebb and flood
        elif len(u_flood_times)==len(u_ebb_times)+1:
            UT1=(u_flood_heights[:-1]+u_ebb_heights)/2
            UT2=(u_flood_heights[1:]+u_ebb_heights)/2
            t_UT1=u_flood_times[:-1]+((u_ebb_times-u_flood_times[:-1])/2)
            t_UT2=u_flood_times[1:]+((u_ebb_times-u_flood_times[1:])/2)
        else:
            print('problem with lengths of flood and ebb peak arrays')
    elif u_flood_times[0]>u_ebb_times[0]:
        if len(u_flood_times)==len(u_ebb_times):
            UT1=(u_flood_heights[:-1]+u_ebb_heights[1:])/2
            UT2=(u_flood_heights+u_ebb_heights)/2
            t_UT1=u_flood_times[:-1]+((u_ebb_times[1:]-u_flood_times[:-1])/2)
            t_UT2=u_flood_times+((u_ebb_times-u_flood_times)/2)
        elif len(u_flood_times)==len(u_ebb_times)-1:
            UT1=(u_flood_heights+u_ebb_heights[1:])/2
            UT2=(u_flood_heights+u_ebb_heights[:-1])/2
            t_UT1=u_flood_times+((u_ebb_times[1:]-u_flood_times)/2)
            t_UT2=u_flood_times+((u_ebb_times[:-1]-u_flood_times)/2)
        else:
            print('problem with lengths of flood and ebb peak arrays')
    #now pick the bigger of UT1 and UT2 and find the corresponding time
    #want to choose only in a range of 1 spring-neap centered around t_spring
    if np.max(UT1[np.abs(t_UT1-t_spring)<np.timedelta64(8,'D')])>np.max(UT2[np.abs(t_UT2-t_spring)<np.timedelta64(8,'D')]): #choose between UT1 and UT2
        t_nearmax=t_UT1[np.abs(t_UT1-t_spring)<np.timedelta64(8,'D')] #select around t_spring
        UT_nearmax=UT1[np.abs(t_UT1-t_spring)<np.timedelta64(8,'D')] #select around t_spring
        t_max=t_nearmax[np.argmax(UT_nearmax)]
        UT_max_v2=UT_nearmax[np.argmax(UT_nearmax)]
    else:
        t_nearmax=t_UT2[np.abs(t_UT2-t_spring)<np.timedelta64(8,'D')] #select around t_spring
        UT_nearmax=UT2[np.abs(t_UT2-t_spring)<np.timedelta64(8,'D')] #select around t_spring
        t_max=t_nearmax[np.argmax(UT_nearmax)]
        UT_max_v2=UT_nearmax[np.argmax(UT_nearmax)]
    t_max_flood = u_flood_times[np.argmin(np.abs(u_flood_times-t_max))]
    t_max_ebb = u_ebb_times[np.argmin(np.abs(u_ebb_times-t_max))]

    #FINDING START OF FLOOD AND EBB (to use for particle tracking releases)
    qnet_pre_spring_flood=qnet_pf[((t_spring_flood-ot_pf)<np.timedelta64(7,'h')) & (ot_pf<=t_spring_flood)] #select qnet for 7h before max flood
    ot_pre_spring_flood=ot_pf[((t_spring_flood-ot_pf)<np.timedelta64(7,'h')) & (ot_pf<=t_spring_flood)] #select times for 7h before max flood
    t_start_spring_flood=ot_pre_spring_flood[np.argmax(qnet_pre_spring_flood>0)] #pick first positive time
    qnet_start_spring_flood=qnet_pre_spring_flood[np.argmax(qnet_pre_spring_flood>0)] #pick corresponding qnet (useful for plotting)

    qnet_pre_spring_ebb=qnet_pf[((t_spring_ebb-ot_pf)<np.timedelta64(7,'h')) & (ot_pf<=t_spring_ebb)]
    ot_pre_spring_ebb=ot_pf[((t_spring_ebb-ot_pf)<np.timedelta64(7,'h')) & (ot_pf<=t_spring_ebb)]
    t_start_spring_ebb=ot_pre_spring_ebb[np.argmax(qnet_pre_spring_ebb<0)] #for ebb pick first negative time
    qnet_start_spring_ebb=qnet_pre_spring_ebb[np.argmax(qnet_pre_spring_ebb<0)]

    qnet_pre_neap_flood=qnet_pf[((t_neap_flood-ot_pf)<np.timedelta64(7,'h')) & (ot_pf<=t_neap_flood)]
    ot_pre_neap_flood=ot_pf[((t_neap_flood-ot_pf)<np.timedelta64(7,'h')) & (ot_pf<=t_neap_flood)]
    t_start_neap_flood=ot_pre_neap_flood[np.argmax(qnet_pre_neap_flood>0)]
    qnet_start_neap_flood=qnet_pre_neap_flood[np.argmax(qnet_pre_neap_flood>0)]

    qnet_pre_neap_ebb=qnet_pf[((t_neap_ebb-ot_pf)<np.timedelta64(7,'h')) & (ot_pf<=t_neap_ebb)]
    ot_pre_neap_ebb=ot_pf[((t_neap_ebb-ot_pf)<np.timedelta64(7,'h')) & (ot_pf<=t_neap_ebb)]
    t_start_neap_ebb=ot_pre_neap_ebb[np.argmax(qnet_pre_neap_ebb<0)]
    qnet_start_neap_ebb=qnet_pre_neap_ebb[np.argmax(qnet_pre_neap_ebb<0)]

    #calculate velocities and excursion
    unet_spring_flood=unet[ot==t_spring_flood]
    unet_spring_ebb=unet[ot==t_spring_ebb]
    unet_neap_flood=unet[ot==t_neap_flood]
    unet_neap_ebb=unet[ot==t_neap_ebb]
    unet_max_flood=unet[ot==t_max_flood]
    unet_max_ebb=unet[ot==t_max_ebb]

    UT_spring=(unet_spring_flood-unet_spring_ebb)/2
    UT_neap=(unet_neap_flood-unet_neap_ebb)/2
    UT_max=(unet_max_flood-unet_max_ebb)/2

    # calculate freshwater Froude number
    UR=QR/A_sect
    c=np.sqrt(beta*g*socn*H_thalweg)
    Frf=UR/c

    # # calculate Utidal
    # qnet=ds2['qnet'] #for making poster
    # qprism=zfun.lowpass(np.abs(ds2['qnet'].values-zfun.lowpass(ds2['qnet'].values, f='godin')), f='godin')/2
    # ds2time=ds2.time
    # unet=ds2['qnet']/A_sect
    # UT_spring=(unet.sel(time=t_spring_flood)-unet.sel(time=t_spring_ebb))/2
    # UT_neap=(unet.sel(time=t_neap_flood)-unet.sel(time=t_neap_ebb))/2
    
    # calculate M
    N0=np.sqrt((beta*g*socn)/H_thalweg)
    M_spring=np.sqrt((Cd*UT_spring*UT_spring)/(omega*N0*H_thalweg*H_thalweg))
    M_neap=np.sqrt((Cd*UT_neap*UT_neap)/(omega*N0*H_thalweg*H_thalweg))

    # calculate tidal excursion
    TE_spring = (UT_spring*Ts)/np.pi
    TE_neap = (UT_neap*Ts)/np.pi
    TE_max = (UT_max*Ts)/np.pi

    # get section longitude
    # out_fn = ext_fn.replace('.nc','.p')
    # one_section_df = pd.read_pickle(c_dir / out_fn)
    # sect_lon.append = one_section_df.loc[0,'x']
    # sect_tick.append(sect_label)


    # #other variables
    # dA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() #shape NT,NZ,NX
    # H = np.sum(ds['DZ'].to_numpy(),axis=1) #shape NT,NX
    # q = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['vel'].to_numpy() #shape NT,NZ,NX
    # sdA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['salt'].to_numpy() #shape NT,NZ,NX
    # # V['q'] = q
    # NT, NZ, NX = q.shape
    # A = np.sum(dA, axis=(1,2)) #shape NT

    #PLOT ON PARAM SPACE
    ax.scatter(M_spring,Frf,s=None,c='tab:green',label='Spring')
    ax.text(M_spring,Frf,sect_label, ha='left',va='top',fontsize=14,c='tab:green')
    ax.scatter(M_neap,Frf,s=None,c='tab:purple',label='Neap')
    ax.text(M_neap,Frf,sect_label, ha='left',va='top',fontsize=14,c='tab:purple')

    if Ldir['testing']:
        if sect_label=='b3':
            ax.legend(loc='upper left')
    else:
        if sect_label=='a1':
            ax.legend(loc='upper left')

    #PLOT TIDE TIMESERIES
    fig1, [ax1,ax2] = plt.subplots(2, 1, sharex=True, figsize=(20,10))
    ax1.plot(ot_pf, qprism_pf)
    ax1.axvline(t_spring,c='tab:green',label='Spring')
    ax1.axvline(t_neap,c='tab:purple',label='Neap')
    ax1.set_title('Qprism')
    ax2.plot(ot_pf,qnet_pf)
    ax2.axvline(t_spring_flood,c='tab:green',linestyle='--',label='Spring flood')
    ax2.axvline(t_spring_ebb,c='tab:green',linestyle=':',label='Spring ebb')
    ax2.axvline(t_neap_flood,c='tab:purple',linestyle='--',label='Neap flood')
    ax2.axvline(t_neap_ebb,c='tab:purple',linestyle=':',label='Neap ebb')
    ax2.axvline(t_max_flood,c='tab:orange',linestyle='--',label='Max UT flood')
    ax2.axvline(t_max_ebb,c='tab:orange',linestyle=':',label='Max UT ebb')
    ax2.scatter(t_start_spring_flood, qnet_start_spring_flood, c='tab:green', marker='s', label='Start of spring flood')
    ax2.scatter(t_start_spring_ebb, qnet_start_spring_ebb, c='tab:green', marker='^', label='Start of spring ebb')
    ax2.scatter(t_start_neap_flood, qnet_start_neap_flood, c='tab:purple', marker='s', label='Start of neap flood')
    ax2.scatter(t_start_neap_ebb, qnet_start_neap_ebb, c='tab:purple', marker='^', label='Start of neap ebb')
    ax2.set_title('qnet')
    ax1.set_xlim(nearday-np.timedelta64(14,'D'),nearday+np.timedelta64(14,'D'))

    ax1.legend()
    ax2.legend()
    ax1.grid(True)
    ax2.grid(True)
    fig1.savefig(out_dir / ('tide_plot.png'))

    ds.close()

    #dict of useful variables
    TE = dict()
    TE['t_spring']=t_spring
    TE['t_neap']=t_neap
    TE['t_max']=t_max
    TE['t_spring_flood']=t_spring_flood
    TE['t_spring_ebb']=t_spring_ebb
    TE['t_neap_flood']=t_neap_flood
    TE['t_neap_ebb']=t_neap_ebb
    TE['t_max_flood']=t_max_flood
    TE['t_max_ebb']=t_max_ebb
    TE['t_start_spring_flood']=t_start_spring_flood
    TE['t_start_spring_ebb']=t_start_spring_ebb
    TE['t_start_neap_flood']=t_start_neap_flood
    TE['t_start_neap_ebb']=t_start_neap_ebb
    TE['unet_spring_flood']=unet_spring_flood
    TE['unet_spring_ebb']=unet_spring_ebb
    TE['unet_neap_flood']=unet_neap_flood
    TE['unet_neap_ebb']=unet_neap_ebb
    TE['unet_max_flood']=unet_max_flood
    TE['unet_max_ebb']=unet_max_ebb
    TE['UT_spring']=UT_spring
    TE['UT_neap']=UT_neap
    TE['UT_max']=UT_max
    TE['UT_max_v2']=UT_max_v2
    TE['TE_spring']=TE_spring
    TE['TE_neap']=TE_neap
    TE['TE_max']=TE_max
    
    TE['UR']=UR
    TE['c']=c
    TE['Frf']=Frf
    TE['N0']=N0
    TE['M_spring']=M_spring
    TE['M_neap']=M_neap
    TE['nearday']=nearday

    out_fn = 'TE_'+sect_label+'.p'
    pickle.dump(TE, open(out_dir/out_fn,'wb'))

    print('t_spring=')
    print(t_spring)
    print('t_neap=')
    print(t_neap)
    print('t_max=')
    print(t_max)
    print('t_spring_flood=')
    print(t_spring_flood)
    print('t_spring_ebb=')
    print(t_spring_ebb)
    print('t_neap_flood=')
    print(t_neap_flood)
    print('t_neap_ebb=')
    print(t_neap_ebb)
    print('UT_spring=')
    print(UT_spring)
    print('UT_neap=')
    print(UT_neap)
    print('UT_max=')
    print(UT_max)
    print('UT_max_v2=')
    print(UT_max_v2)
    print('TE_spring (km)=')
    print(TE_spring/1000)
    print('TE_neap (km)=')
    print(TE_neap/1000)
    print('TE_max (km)=')
    print(TE_max/1000)


    print('  elapsed time for section = %d seconds' % (time()-tt0))
    sys.stdout.flush()

#axes limits to match G&M 2014 plot
#remove axes limits temporarily
# ax.set_xlim(1.5,2)
# ax.set_ylim(1e-4,1e0)

ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'$M$')
ax.set_ylabel(r'$Fr_f$')
ax.grid(True)

fig.suptitle('Idealized model sections in estuarine parameter space')
fig.savefig(out_dir / ('param_plot.png'))
plt.close()

# #Plot tidal excursion
# TE_spring=np.asarray(TE_spring)
# TE_neap=np.asarray(TE_neap)
# #sect_lon=np.asarray(sect_lon)

# pfun.start_plot(fs=fs, figsize=(15,15))
# fig, ax = plt.subplots(1, 1)

# xplot=np.arange(0,len(TE_spring))
# ax.plot(xplot,TE_spring/1000,'-o',c='tab:green',label='Spring')
# ax.plot(xplot,TE_neap/1000,'-o',c='tab:purple',label='Neap')
# ax.set_ylim(bottom=0)
# ax.set_xlim(np.min(xplot),np.max(xplot))
# ax.set_xticks(xplot,sect_tick)

# #ax.set_xlim(1.5,2)
# #ax.set_ylim(1e-4,1e0)

# ax.grid(True)
# ax.set_xlabel('Section')
# ax.set_ylabel('Tidal excursion [km]')
# ax.legend(loc='upper right')

# fig.suptitle('Tidal excursion along estuary')
# plt.savefig(out_dir / ('excur_plot.png'))
# plt.close()


print('\nTotal elapsed time for calculation and plotting = %d seconds' % (time()-tt00))



