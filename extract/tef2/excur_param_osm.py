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
import matplotlib.patches as patches
from cmocean import cm

from lo_tools import plotting_functions as pfun

from lo_tools import Lfun, zrfun, zfun
from lo_tools import extract_argfun as exfun
# Ldir = exfun.intro() # this handles the argument passing
# use the arg -sect_name to pick section for plots

# import tef_fun


# gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
# tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

# sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
# sect_df = pd.read_pickle(sect_df_fn)
# c_dir = Ldir['LOo'] / 'extract' / 'tef2' / ('sections' + '_' + gctag)

# out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
# in_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
# in_dir2 = out_dir0 / ('processed_' + Ldir['ds0'] + '_' + Ldir['ds1']) #in_dir for processed to get qnet
# out_dir = out_dir0 / ('excur_param_plots_' + Ldir['ds0'] + '_' + Ldir['ds1'])
# Lfun.make_dir(out_dir, clean=True)

# sect_list = [item.name for item in in_dir.glob('*.nc')]
# if Ldir['testing']:
#     sect_list = ['b1.nc','b2.nc','b3.nc','b4.nc','b5.nc']
# #sect_list = [Ldir['sect_name']+'.nc'] #section to plot #not sure if this syntax is correct
    
# # make vn_list by inspecting the first section
# ds = xr.open_dataset(in_dir / sect_list[0])
# vn_list = [item for item in ds.data_vars \
#     if (len(ds[item].dims) == 3) and (item not in ['vel','DZ'])]
# ds.close()


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

#values for current model
UT_spring=0.5
UT_neap=0.25
H_sill=50
H_basin=200
W_sill=4000
W_basin=8000
QR=1000
#calculate freshwater Froude number
A_sect = H_sill*W_sill
UR=QR/A_sect
c=np.sqrt(beta*g*socn*H_sill)
Frf=UR/c
#calculate M
N0=np.sqrt((beta*g*socn)/H_sill)
M_spring=np.sqrt((Cd*UT_spring*UT_spring)/(omega*N0*H_sill*H_sill))
M_neap=np.sqrt((Cd*UT_neap*UT_neap)/(omega*N0*H_sill*H_sill))
#calculate tidal excursion
TE_spring=(UT_spring*Ts)/np.pi
TE_neap=(UT_neap*Ts)/np.pi
print('\nspring tidal excursion:')
print(TE_spring/1000)
print('\nneap tidal excursion')
print(TE_neap/1000)
#labels
current_label=' Current'

# #vary QR
# QRvary=np.array([1000,2000,3000,4000])
# UR_QRvary=QRvary/A_sect
# Frf_QRvary=UR_QRvary/c
# M_spring_QRvary=M_spring*np.ones(len(QRvary))
# M_neap_QRvary=M_neap*np.ones(len(QRvary))
# #labels
# QR_labels=['',' 2x QR',' 3x QR',' 4x QR']

# #vary tidal amplitude
# TA_factor=np.array([1,2,3])
# UT_spring_TAvary=UT_spring*TA_factor
# UT_neap_TAvary=UT_neap*TA_factor
# M_spring_TAvary=np.sqrt((Cd*UT_spring_TAvary*UT_spring_TAvary)/(omega*N0*H_sill*H_sill))
# M_neap_TAvary=np.sqrt((Cd*UT_neap_TAvary*UT_neap_TAvary)/(omega*N0*H_sill*H_sill))
# Frf_TAvary=Frf*np.ones(len(TA_factor))
# #labels
# TA_labels=['',' 2x tide',' 3x tide']

# #vary H
# H_factor=np.array([1, 1/2, 1/3, 1/4])
# Hvary=H_sill*H_factor
# UT_spring_Hvary=UT_spring*(1/H_factor)
# UT_neap_Hvary=UT_neap*(1/H_factor)
# N0_Hvary=np.sqrt((beta*g*socn)/Hvary)
# M_spring_Hvary=np.sqrt((Cd*UT_spring_Hvary*UT_spring_Hvary)/(omega*N0_Hvary*Hvary*Hvary))
# M_neap_Hvary=np.sqrt((Cd*UT_neap_Hvary*UT_neap_Hvary)/(omega*N0_Hvary*Hvary*Hvary))
# A_sect_Hvary=Hvary*W_sill
# c_Hvary=np.sqrt(beta*g*socn*Hvary)
# UR_Hvary=QR/A_sect_Hvary
# Frf_Hvary=UR_Hvary/c_Hvary
# #labels
# H_labels=['',' 1/2x H',' 1/3x H',' 1/4x H']

# #vary W
# W_factor=np.array([1, 1/2, 1/3, 1/4])
# Wvary=W_sill*W_factor
# UT_spring_Wvary=UT_spring*(1/W_factor)
# UT_neap_Wvary=UT_neap*(1/W_factor)
# M_spring_Wvary=np.sqrt((Cd*UT_spring_Wvary*UT_spring_Wvary)/(omega*N0*H_sill*H_sill))
# M_neap_Wvary=np.sqrt((Cd*UT_neap_Wvary*UT_neap_Wvary)/(omega*N0*H_sill*H_sill))
# A_sect_Wvary=H_sill*Wvary
# UR_Wvary=QR/A_sect_Wvary
# Frf_Wvary=UR_Wvary/c
# #labels
# W_labels=['',' 1/2x W',' 1/3x W',' 1/4x W']


#start parameter space plot
plt.close('all')
fs = 14
pfun.start_plot(fs=fs, figsize=(15,15))
fig, ax = plt.subplots(1, 1)

# #plot current model on param space
# ax.scatter(M_spring,Frf,s=None,c='tab:green',label='Spring')
# ax.text(M_spring,Frf,current_label, ha='right',va='top',fontsize=14,c='tab:green')
# ax.scatter(M_neap,Frf,s=None,c='tab:purple',label='Neap')
# ax.text(M_neap,Frf,current_label, ha='right',va='top',fontsize=14,c='tab:purple')
# ax.legend(loc='upper left')





# #plot Q_vary lines
# ax.plot(M_spring_QRvary,Frf_QRvary,'-.o',c='tab:green',label='Neap')
# ax.plot(M_neap_QRvary,Frf_QRvary,'-.o',c='tab:purple',label='Spring')
# for i in range(len(QRvary)):
#     ax.text(M_spring_QRvary[i],Frf_QRvary[i],QR_labels[i], ha='right',va='top',fontsize=14,c='tab:green')
#     ax.text(M_neap_QRvary[i],Frf_QRvary[i],QR_labels[i], ha='right',va='top',fontsize=14,c='tab:purple')

# #plot TA_vary lines
# ax.plot(M_spring_TAvary,Frf_TAvary,':o',c='tab:green',label='Neap')
# ax.plot(M_neap_TAvary,Frf_TAvary,':o',c='tab:purple',label='Spring')
# for i in range(len(TA_factor)):
#     ax.text(M_spring_TAvary[i],Frf_TAvary[i],TA_labels[i], ha='right',va='bottom',fontsize=14,c='tab:green')
#     ax.text(M_neap_TAvary[i],Frf_TAvary[i],TA_labels[i], ha='right',va='bottom',fontsize=14,c='tab:purple')

# #plot H_vary lines
# ax.plot(M_spring_Hvary,Frf_Hvary,'-o',c='tab:green',label='Neap')
# ax.plot(M_neap_Hvary,Frf_Hvary,'-o',c='tab:purple',label='Spring')
# for i in range(len(H_factor)):
#     ax.text(M_spring_Hvary[i],Frf_Hvary[i],H_labels[i], ha='left',va='top',fontsize=14,c='tab:green')
#     ax.text(M_neap_Hvary[i],Frf_Hvary[i],H_labels[i], ha='left',va='top',fontsize=14,c='tab:purple')

# #plot W_vary lines
# ax.plot(M_spring_Wvary,Frf_Wvary,'--o',c='tab:green',label='Neap')
# ax.plot(M_neap_Wvary,Frf_Wvary,'--o',c='tab:purple',label='Spring')
# for i in range(len(W_factor)):
#     ax.text(M_spring_Wvary[i],Frf_Wvary[i],W_labels[i], ha='right',va='bottom',fontsize=14,c='tab:green')
#     ax.text(M_neap_Wvary[i],Frf_Wvary[i],W_labels[i], ha='right',va='bottom',fontsize=14,c='tab:purple')

# Known estuaries
PS_bottom=1.18e-4
PS_top=3.16e-4
PS_left=0.15
PS_right=0.325
PS_height=PS_top-PS_bottom
PS_width=PS_right-PS_left
PS_rect=plt.Rectangle((PS_left,PS_bottom), PS_width, PS_height,fill=False , edgecolor='tab:blue')
ax.add_patch(PS_rect)
ax.text(PS_left,PS_top,'Puget Sound', ha='left',va='top',fontsize=14,c='tab:blue')

CB_bottom=1.54e-3
CB_top=1.33e-2
CB_left=0.308
CB_right=0.548
CB_height=CB_top-CB_bottom
CB_width=CB_right-CB_left
CB_rect=plt.Rectangle((CB_left,CB_bottom), CB_width, CB_height,fill=False , edgecolor='tab:olive')
ax.add_patch(CB_rect)
ax.text(CB_left,CB_top,'Chesapeake Bay', ha='left',va='top',fontsize=14,c='tab:olive')

HR_bottom=6.49e-3
HR_top=1.1e-1
HR_left=0.53
HR_right=0.75
HR_height=HR_top-HR_bottom
HR_width=HR_right-HR_left
HR_rect=plt.Rectangle((HR_left,HR_bottom), HR_width, HR_height,fill=False , edgecolor='tab:brown')
ax.add_patch(HR_rect)
ax.text(HR_left,HR_top,'Hudson River', ha='left',va='top',fontsize=14,c='tab:brown')

#Admiralty inlet
AI_UTspring=2
AI_UTneap=1
AI_H=65
AI_socn=32
AI_QRlow=800
AI_QRhigh=1500
AI_W=6000
AI_Ls=30000
#calculate freshwater Froude number
AI_A_sect = AI_H*AI_W
AI_URlow=AI_QRlow/AI_A_sect
AI_URhigh=AI_QRhigh/AI_A_sect
AI_c=np.sqrt(beta*g*AI_socn*AI_H)
AI_Frflow=AI_URlow/AI_c
AI_Frfhigh=AI_URhigh/AI_c
#calculate M
AI_N0=np.sqrt((beta*g*AI_socn)/AI_H)
AI_Mspring=np.sqrt((Cd*AI_UTspring*AI_UTspring)/(omega*AI_N0*AI_H*AI_H))
AI_Mneap=np.sqrt((Cd*AI_UTneap*AI_UTneap)/(omega*AI_N0*AI_H*AI_H))
#calculate tidal excursion
AI_TEspring=(AI_UTspring*Ts)/np.pi
AI_TEneap=(AI_UTneap*Ts)/np.pi
print('\nAI spring tidal excursion:')
print(AI_TEspring/1000)
print('\nAI neap tidal excursion')
print(AI_TEneap/1000)
#add to plot
AI_rect=plt.Rectangle((AI_Mneap,AI_Frflow), AI_Mspring-AI_Mneap, AI_Frfhigh-AI_Frflow,fill=False , edgecolor='tab:cyan')
ax.add_patch(AI_rect)
ax.text(AI_Mneap,AI_Frfhigh,'Admiralty Inlet', ha='left',va='top',fontsize=14,c='tab:cyan')

#plot 20km model on param space
M_spring_20=0.545
M_neap_20=0.267
Frf_20=0.00189

ax.scatter(M_spring_20,Frf_20,s=None,c='tab:green',label='Spring')
#ax.text(M_spring_20,Frf_20,'Spring', ha='left',va='center',fontsize=14,c='k')
ax.scatter(M_neap_20,Frf_20,s=None,c='tab:green',label='Neap')
#ax.text(M_neap_20,Frf_20,'Neap', ha='right',va='center',fontsize=14,c='k')
ax.plot([M_neap_20,M_spring_20],[Frf_20,Frf_20],':',c='tab:green')
ax.text((M_spring_20+M_neap_20)/2,Frf_20,'20km sill', ha='center',va='bottom',fontsize=14,c='tab:green')

#plot 5km model on param space
M_spring_5=0.454
M_neap_5=0.217
Frf_5=0.00189

ax.scatter(M_spring_5,Frf_5,s=None,c='r',label='Spring')
# ax.text(M_spring_20,Frf_20,'Spring', ha='left',va='center',fontsize=14,c='k')
ax.scatter(M_neap_5,Frf_5,s=None,c='r',label='Neap')
# ax.text(M_neap_20,Frf_20,'Neap', ha='right',va='center',fontsize=14,c='k')
ax.plot([M_neap_5,M_spring_5],[Frf_5,Frf_5],':',c='r')
ax.text((M_spring_5+M_neap_5)/2,Frf_5,'5km sill', ha='center',va='bottom',fontsize=14,c='r')

#plot 5km model on param space
M_spring_80=1.06
M_neap_80=0.602
Frf_80=0.00189

ax.scatter(M_spring_80,Frf_80,s=None,c='tab:purple',label='Spring')
# ax.text(M_spring_20,Frf_20,'Spring', ha='left',va='center',fontsize=14,c='k')
ax.scatter(M_neap_80,Frf_80,s=None,c='tab:purple',label='Neap')
# ax.text(M_neap_20,Frf_20,'Neap', ha='right',va='center',fontsize=14,c='k')
ax.plot([M_neap_80,M_spring_80],[Frf_80,Frf_80],':',c='tab:purple')
ax.text((M_spring_80+M_neap_80)/2,Frf_80,'80km sill', ha='center',va='bottom',fontsize=14,c='tab:purple')

#format plot
#axes limits to match G&M 2014 plot
#remove axes limits temporarily
# ax.set_xlim(1.5,2)
# ax.set_ylim(1e-4,1e0)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel(r'Mixing parameter $M$')
ax.set_ylabel(r'Freshwater Froude number $Fr_f$')
ax.grid(True)
fig.suptitle('Idealized model in estuarine parameter space')
# plt.savefig(out_dir / ('param_plot.png'))
# plt.close()
plt.show()


# #list for tidal excursion
# TE_spring=[]
# TE_neap=[]
# sect_tick=[]

# for ext_fn in sect_list:
#     tt0 = time()
#     print(ext_fn)
#     sect_label = ext_fn.replace('.nc','')

#     # load extraction fields
#     ds = xr.open_dataset(in_dir / ext_fn)

#     # load processed fields
#     ds2 = xr.open_dataset(in_dir2 / ext_fn)
#     #ax9.plot(ds2['time'],ds2['qnet'])


#     # get variables from section extraction
#     ot = ds['time'].to_numpy() #shape NT
#     zeta = ds['zeta'].to_numpy() #shape NT,NX
#     h = ds['h'].to_numpy() #shape NX
#     dd = ds['dd'].to_numpy() #shape NX
#     dz = ds['DZ'].to_numpy() #shape NT,NZ,NX
#     u = ds['vel'].to_numpy() #shape NT,NZ,NX
#     s = ds['salt'].to_numpy() #shape NT,NZ,NX

#     # calculate section variables
#     H_thalweg=np.max(h) #thalweg depth (use as depth for N0 and M)
#     A_sect=np.sum(h*dd) #section area (not time varying)

#     # calculate freshwater Froude number
#     UR=QR/A_sect
#     c=np.sqrt(beta*g*socn*H_thalweg)
#     Frf=UR/c

#     # calculate Utidal
#     unet=ds2['qnet']/A_sect
#     UT_spring=(unet.sel(time=t_spring_flood)-unet.sel(time=t_spring_ebb))/2
#     UT_neap=(unet.sel(time=t_neap_flood)-unet.sel(time=t_neap_ebb))/2
    
#     # calculate M
#     N0=np.sqrt((beta*g*socn)/H_thalweg)
#     M_spring=np.sqrt((Cd*UT_spring*UT_spring)/(omega*N0*H_thalweg*H_thalweg))
#     M_neap=np.sqrt((Cd*UT_neap*UT_neap)/(omega*N0*H_thalweg*H_thalweg))

#     # calculate tidal excursion
#     TE_spring.append((UT_spring*Ts)/np.pi)
#     TE_neap.append((UT_neap*Ts)/np.pi)

#     # get section longitude
#     # out_fn = ext_fn.replace('.nc','.p')
#     # one_section_df = pd.read_pickle(c_dir / out_fn)
#     # sect_lon.append = one_section_df.loc[0,'x']
#     sect_tick.append(sect_label)


#     # #other variables
#     # dA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() #shape NT,NZ,NX
#     # H = np.sum(ds['DZ'].to_numpy(),axis=1) #shape NT,NX
#     # q = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['vel'].to_numpy() #shape NT,NZ,NX
#     # sdA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['salt'].to_numpy() #shape NT,NZ,NX
#     # # V['q'] = q
#     # NT, NZ, NX = q.shape
#     # A = np.sum(dA, axis=(1,2)) #shape NT

#     #PLOT ON PARAM SPACE
#     ax.scatter(M_spring,Frf,s=None,c='tab:green',label='Neap')
#     ax.text(M_spring,Frf,sect_label, ha='left',va='top',fontsize=14,c='tab:green')
#     ax.scatter(M_neap,Frf,s=None,c='tab:purple',label='Spring')
#     ax.text(M_neap,Frf,sect_label, ha='left',va='top',fontsize=14,c='tab:purple')

#     if Ldir['testing']:
#         if sect_label=='b1':
#             ax.legend(loc='upper left')
#     else:
#         if sect_label=='a1':
#             ax.legend(loc='upper left')

#     ds.close()
#     ds2.close()

#     print('  elapsed time for section = %d seconds' % (time()-tt0))
#     sys.stdout.flush()



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

# print('\nTotal elapsed time for calculation and plotting = %d seconds' % (time()-tt00))



