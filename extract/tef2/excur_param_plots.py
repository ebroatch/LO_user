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
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from cmocean import cm

from lo_tools import plotting_functions as pfun

from lo_tools import Lfun, zrfun, zfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing
# use the arg -sect_name to pick section for plots

# import tef_fun


gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)
c_dir = Ldir['LOo'] / 'extract' / 'tef2' / ('sections' + '_' + gctag)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
in_dir2 = out_dir0 / ('processed_' + Ldir['ds0'] + '_' + Ldir['ds1']) #in_dir for processed to get qnet
out_dir = out_dir0 / ('excur_param_plots_' + Ldir['ds0'] + '_' + Ldir['ds1'])
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
# t_spring_ebb = pd.Timestamp('2020-07-01 04:00:00') #old sill3 model
# t_spring_flood = pd.Timestamp('2020-07-01 10:00:00')
# t_neap_ebb = pd.Timestamp('2020-07-08 10:00:00')
# t_neap_flood = pd.Timestamp('2020-07-08 16:00:00')

# t_spring = pd.Timestamp('2020-07-01 07:00:00')
# t_neap = pd.Timestamp('2020-07-08 14:00:00')

#t_spring_ebb = pd.Timestamp('2020-10-12 17:00:00') #20km model
t_spring_ebb = pd.Timestamp('2020-10-12 16:00:00') #5km model
t_spring_flood = pd.Timestamp('2020-10-12 10:00:00')
t_neap_ebb = pd.Timestamp('2020-10-05 11:00:00')
t_neap_flood = pd.Timestamp('2020-10-05 04:00:00')

#list for tidal excursion
TE_spring=[]
TE_neap=[]
sect_tick=[]

#start parameter space plot
plt.close('all')
fs = 14
pfun.start_plot(fs=fs, figsize=(15,15))
fig, ax = plt.subplots(1, 1)

for ext_fn in sect_list:
    tt0 = time()
    print(ext_fn)
    sect_label = ext_fn.replace('.nc','')

    # load extraction fields
    ds = xr.open_dataset(in_dir / ext_fn)

    # load processed fields
    ds2 = xr.open_dataset(in_dir2 / ext_fn)
    #ax9.plot(ds2['time'],ds2['qnet'])


    # get variables from section extraction
    ot = ds['time'].to_numpy() #shape NT
    zeta = ds['zeta'].to_numpy() #shape NT,NX
    h = ds['h'].to_numpy() #shape NX
    dd = ds['dd'].to_numpy() #shape NX
    dz = ds['DZ'].to_numpy() #shape NT,NZ,NX
    u = ds['vel'].to_numpy() #shape NT,NZ,NX
    s = ds['salt'].to_numpy() #shape NT,NZ,NX

    # calculate section variables
    H_thalweg=np.max(h) #thalweg depth (use as depth for N0 and M)
    A_sect=np.sum(h*dd) #section area (not time varying)

    # calculate freshwater Froude number
    UR=QR/A_sect
    c=np.sqrt(beta*g*socn*H_thalweg)
    Frf=UR/c

    # calculate Utidal
    qnet=ds2['qnet'] #for making poster
    qprism=zfun.lowpass(np.abs(ds2['qnet'].values-zfun.lowpass(ds2['qnet'].values, f='godin')), f='godin')/2
    ds2time=ds2.time
    unet=ds2['qnet']/A_sect
    UT_spring=(unet.sel(time=t_spring_flood)-unet.sel(time=t_spring_ebb))/2
    UT_neap=(unet.sel(time=t_neap_flood)-unet.sel(time=t_neap_ebb))/2
    
    # calculate M
    N0=np.sqrt((beta*g*socn)/H_thalweg)
    M_spring=np.sqrt((Cd*UT_spring*UT_spring)/(omega*N0*H_thalweg*H_thalweg))
    M_neap=np.sqrt((Cd*UT_neap*UT_neap)/(omega*N0*H_thalweg*H_thalweg))

    # calculate tidal excursion
    TE_spring.append((UT_spring*Ts)/np.pi)
    TE_neap.append((UT_neap*Ts)/np.pi)

    # get section longitude
    # out_fn = ext_fn.replace('.nc','.p')
    # one_section_df = pd.read_pickle(c_dir / out_fn)
    # sect_lon.append = one_section_df.loc[0,'x']
    sect_tick.append(sect_label)


    # #other variables
    # dA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() #shape NT,NZ,NX
    # H = np.sum(ds['DZ'].to_numpy(),axis=1) #shape NT,NX
    # q = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['vel'].to_numpy() #shape NT,NZ,NX
    # sdA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['salt'].to_numpy() #shape NT,NZ,NX
    # # V['q'] = q
    # NT, NZ, NX = q.shape
    # A = np.sum(dA, axis=(1,2)) #shape NT

    #PLOT ON PARAM SPACE
    ax.scatter(M_spring,Frf,s=None,c='tab:green',label='Neap')
    ax.text(M_spring,Frf,sect_label, ha='left',va='top',fontsize=14,c='tab:green')
    ax.scatter(M_neap,Frf,s=None,c='tab:purple',label='Spring')
    ax.text(M_neap,Frf,sect_label, ha='left',va='top',fontsize=14,c='tab:purple')

    if Ldir['testing']:
        if sect_label=='b1':
            ax.legend(loc='upper left')
    else:
        if sect_label=='a1':
            ax.legend(loc='upper left')

    ds.close()
    ds2.close()

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
plt.savefig(out_dir / ('param_plot.png'))
plt.close()

#Plot tidal excursion
TE_spring=np.asarray(TE_spring)
TE_neap=np.asarray(TE_neap)
#sect_lon=np.asarray(sect_lon)

pfun.start_plot(fs=fs, figsize=(15,15))
fig, ax = plt.subplots(1, 1)

xplot=np.arange(0,len(TE_spring))
ax.plot(xplot,TE_spring/1000,'-o',c='tab:green',label='Spring')
ax.plot(xplot,TE_neap/1000,'-o',c='tab:purple',label='Neap')
ax.set_ylim(bottom=0)
ax.set_xlim(np.min(xplot),np.max(xplot))
ax.set_xticks(xplot,sect_tick)

#ax.set_xlim(1.5,2)
#ax.set_ylim(1e-4,1e0)

ax.grid(True)
ax.set_xlabel('Section')
ax.set_ylabel('Tidal excursion [km]')
ax.legend(loc='upper right')

fig.suptitle('Tidal excursion along estuary')
plt.savefig(out_dir / ('excur_plot.png'))
plt.close()

print('\nTotal elapsed time for calculation and plotting = %d seconds' % (time()-tt00))



