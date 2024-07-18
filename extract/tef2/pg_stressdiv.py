"""
Code to make s(z) on multiple sections.

run dsdx_spring_neap -gtx cas6_v0_live -ctag c0 -0 2018.01.01 -1 2018.12.31

"""

import sys
import xarray as xr
import numpy as np
import pickle
from time import time
import pandas as pd
from scipy.stats import binned_statistic
import seawater as sw

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from lo_tools import plotting_functions as pfun
import matplotlib.colors as colors
from cmocean import cm

from lo_tools import Lfun, zrfun, zfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

# output location
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
out_dir = out_dir0 / ('dsdx_plots_' + Ldir['ds0'] + '_' + Ldir['ds1'])
#out_dir = Ldir['parent'] / 'LPM_output' / 'extract'/ 'tef_exdyn'
Lfun.make_dir(out_dir)

# gctag and location of tef2 section definitions
gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

# get sect_df with the section point locations
sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

# get tide info from the tide excursion calculator
excur_dir = out_dir0 / ('tide_excursion_' + Ldir['ds0'] + '_' + Ldir['ds1'])
te_fn = excur_dir / ('TE_b3.p') #could change if using other sections
TE = pd.read_pickle(te_fn)

# get the grid file
gds = xr.open_dataset(Ldir['grid'] / 'grid.nc')
lou = gds.lon_u[0,:].values
lau = gds.lat_u[:,0].values
lov = gds.lon_v[0,:].values
lav = gds.lat_v[:,0].values
lor = gds.lon_rho.values
lar = gds.lat_rho.values
plon, plat = pfun.get_plon_plat(lor,lar)
hh = gds.h.values
maskr = gds.mask_rho.values
zm = -np.ma.masked_where(maskr==0, hh)

# create the dict S
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)

# where to find the extracted sections
in_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = in_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
in_dir2 = in_dir0 / ('extractions_stress_' + Ldir['ds0'] + '_' + Ldir['ds1'])
in_dir3 = in_dir0 / ('extractions_rpm_' + Ldir['ds0'] + '_' + Ldir['ds1'])

# define the list of sections to work on
sect_list = [item.name for item in in_dir.glob('*.nc')]
sect_list = [item.replace('.nc','') for item in sect_list]

# Define sections to work on.
# Generally choose [seaward, landward]
#sect_list = ['b3'] #TESTING
sect_list = ['b1','b2','b3','b4','b5'] #SHORT LIST SILL ONLY
#sect_list = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5'] #FULL LIST
#sect_list = ['ai1','ai2','ai4','ai5','ai6','ai7'] # AI North to South
#sect_list = ['sog7','sog6','sog5','sog4','sog3','sog2'] # SoG North to South
#sect_list = ['jdf1','jdf2','jdf3','jdf4','sji1','sji4'] # JdF to Haro Strait
color_list1 = [plt.cm.tab20(6),plt.cm.tab20(2),plt.cm.tab20(4),plt.cm.tab20(18),plt.cm.tab20(0)]
color_list2 = [plt.cm.tab20(7),plt.cm.tab20(3),plt.cm.tab20(5),plt.cm.tab20(19),plt.cm.tab20(1)]

# make vn_list by inspecting the first section
ds = xr.open_dataset(in_dir / (sect_list[0] + '.nc'))
vn_list = [item for item in ds.data_vars \
    if (len(ds[item].dims) == 3) and (item not in ['vel','DZ'])]
ds.close()

print('\nCalculating s(z) on mutiple sections')
print(str(in_dir))

tt00 = time()

stz_dict = dict()
dustrdztz_dict = dict()
zeta_dict = dict()
lon_dict = dict()
lat_dict = dict()
lon_vec_dict = dict()
lat_vec_dict = dict()
H_dict = dict()
Qprism_dict = dict()

pg_center_neap_dict=dict()
pg_center_spring_dict=dict()
dustrdz_center_neap_dict=dict()
dustrdz_center_spring_dict=dict()
zr_center_dict=dict()

for sn in sect_list:
    tt0 = time()
    print(sn)

    # load fields
    ds = xr.open_dataset(in_dir / (sn + '.nc'))
    salt_hourly = ds.salt.values
    u_hourly = ds.vel.values
    pad = 36
    salt = zfun.lowpass(salt_hourly, f='godin')[pad:-pad+1:24, :]
    NT, N, P = salt.shape

    # # make Qprism
    q = ds['dd'].values * ds['DZ'].values * ds['vel'].values
    qnet = np.nan * np.ones(q.shape[0])
    for tt in range(q.shape[0]):
        qnet[tt] = q[tt,:,:].squeeze().sum()
    qabs = np.abs(qnet)
    qabs_lp = zfun.lowpass(qabs, f='godin')[pad:-pad+1:24]
    Qprism = qabs_lp/2
    Qprism_dict[sn] = Qprism
    
    sdf = sect_df.loc[sect_df.sn==sn,:]
    
    h = ds.h.values
    zr, zw = zrfun.get_z(h, 0*h, S) # packed (z,p)
    zf = zr.flatten() # Note: this does not change with time
    zwf = zw.flatten()
    # what if there is only one p?
    H_dict[sn] = h.mean()
    
    # get section area (ssh=0)
    dz = np.diff(zw,axis=0)
    dzr = np.diff(zr,axis=0)
    A = np.sum(ds['dd'].values * dz)
    
    # Find mean lat and lon (more work than it should be!). #I THINK WE CAN SKIP/CHANGE THIS
    lon_vec = np.concatenate((lou[sdf.loc[(sdf.uv=='u'),'i']],lov[sdf.loc[(sdf.uv=='v'),'i']]))
    lat_vec = np.concatenate((lau[sdf.loc[(sdf.uv=='u'),'j']],lav[sdf.loc[(sdf.uv=='v'),'j']]))
    lon_vec_dict[sn] = lon_vec
    lat_vec_dict[sn] = lat_vec
    lo = np.mean(lon_vec)
    la = np.mean(lat_vec)
    lon_dict[sn] = lo
    lat_dict[sn] = la
    
    # # Then we want to form a time series of s(z)
    # NZ = 100 #ADD MORE BINS???
    # #NZ = 200 #ADD MORE BINS???
    # s_vs_z = np.nan * np.ones((NT,NZ))
    # for tt in range(NT):
    #     sf = salt[tt,:,:].squeeze().flatten()
    #     # scipy.stats.binned_statistic(x, values, statistic='mean', bins=10, range=None)
    #     ##bs = binned_statistic(zf, sf, statistic='mean', bins=NZ, range=(-500,0))
    #     bs = binned_statistic(zf, sf, statistic='mean', bins=NZ, range=(-200,0)) #change to -200 max depth in estuary
    #     s_vs_z[tt,:] = bs.statistic
    # bin_edges = bs.bin_edges

    # load stress fields
    ds2 = xr.open_dataset(in_dir2 / (sn + '.nc'))
    AKv_hourly = ds2.AKv.values
    bustr_hourly = ds2.bustr.values
    RHO0=1023.7 #mean density from BLANK.in file for model
    DZR = ds2['DZR'].values
    pad = 36

    # calculate tidally averaged stress
    dudz_hourly = np.diff(u_hourly,axis=1)/DZR
    ustr_hourly = AKv_hourly*dudz_hourly
    ustr_hourly_full = np.concatenate((bustr_hourly[:,np.newaxis,:]/RHO0,ustr_hourly,np.zeros((bustr_hourly.shape[0],1,bustr_hourly.shape[1]))),axis=1)
    ustr = zfun.lowpass(ustr_hourly_full, f='godin')[pad:-pad+1:24, :] #tidally average and subsample daily
    dustrdz = np.diff(ustr,axis=1)/dz

    # bin stress into vertical bins same as salt
    # ustr_vs_z = np.nan * np.ones((NT,NZ))
    # for tt in range(NT):
    #     ustrf = ustr[tt,:,:].squeeze().flatten()
    #     # scipy.stats.binned_statistic(x, values, statistic='mean', bins=10, range=None)
    #     ##bs = binned_statistic(zf, sf, statistic='mean', bins=NZ, range=(-500,0))
    #     bs2 = binned_statistic(zwf, ustrf, statistic='mean', bins=NZ, range=(-200,0)) #change to -200 max depth in estuary
    #     ustr_vs_z[tt,:] = bs2.statistic
    # bin_edges2 = bs2.bin_edges

    # # bin stress derivative into vertical bins #no more binning
    # dustrdz_vs_z = np.nan * np.ones((NT,NZ))
    # for tt in range(NT):
    #     dustrdzf = dustrdz[tt,:,:].squeeze().flatten()
    #     # scipy.stats.binned_statistic(x, values, statistic='mean', bins=10, range=None)
    #     ##bs = binned_statistic(zf, sf, statistic='mean', bins=NZ, range=(-500,0))
    #     bs2 = binned_statistic(zf, dustrdzf, statistic='mean', bins=NZ, range=(-200,0)) #change to -200 max depth in estuary
    #     dustrdz_vs_z[tt,:] = bs2.statistic
    # bin_edges2 = bs2.bin_edges

    # load rpm fields
    ds3 = xr.open_dataset(in_dir3 / (sn + '.nc'))
    zetarp = ds3.zetarp.values
    zetarm = ds3.zetarm.values
    saltrp = ds3.saltrp.values
    saltrm = ds3.saltrm.values
    DZ = ds3['DZ'].values #note DZ is time varying and dz is not (calculated using zeta=0)
    pad = 36

    # get distance array dx to use with the rpm fields (should be basically the same for each cell)
    lorp = lor[sdf.jrp,sdf.irp]
    lorm = lor[sdf.jrm,sdf.irm]
    larp = lar[sdf.jrp,sdf.irp]
    larm = lar[sdf.jrm,sdf.irm]
    dxrpm = []
    for i in range(len(lorp)):
        dx, ang = sw.dist([larm[i],larp[i]],[lorm[i],lorp[i]],units='km')
        dxrpm.append(dx[0]*1000) #put into m

    # calculate tidally averaged ds/dx and dzeta/dx
    dsdx=zfun.lowpass((saltrp-saltrm)/dxrpm, f='godin')[pad:-pad+1:24, :]
    dzetadx=zfun.lowpass((zetarp-zetarm)/dxrpm, f='godin')[pad:-pad+1:24, :]

    #calculate pressure gradient
    g=9.81
    beta=7.7e-4
    pgzeta = g*dzetadx[:,np.newaxis,:]
    pgs = beta*g*np.flip(np.cumsum(np.flip(dz*dsdx,axis=1),axis=1),axis=1) #need to flip since we want to sum starting from the top which is the last element
    pgscorr = beta*g*( np.flip(np.cumsum(np.flip(dz*dsdx,axis=1),axis=1),axis=1) - 0.5*dz*dsdx)
    pg=pgzeta+pgs 
    pgcorr=pgzeta+pgscorr

    # calculate average zeta across the section
    zeta_avg = np.sum(ds['dd'].values*ds['zeta'].values,axis=1)/np.sum(ds['dd'].values)
    
    ot = ds['time'].to_numpy()
    # do a little massaging of ot
    dti = pd.to_datetime(ot) # a pandas DatetimeIndex with dtype='datetime64[ns]'
    dt = dti.to_pydatetime() # an array of datetimes
    ot = np.array([Lfun.datetime_to_modtime(item) for item in dt])
    ot = ot[pad:-pad+1:24]
    # also make an array of datetimes to save as the ot variable
    otdt = np.array([Lfun.modtime_to_datetime(item) for item in ot])

    #select dustrdz and pg fields closest to spring and neap
    it_neap = zfun.find_nearest_ind(dti[pad:-pad+1:24],TE['t_neap']) #daily subsampling
    it_spring = zfun.find_nearest_ind(dti[pad:-pad+1:24],TE['t_spring'])
    pg_neap=pgcorr[it_neap,:,:]
    pg_spring=pgcorr[it_spring,:,:]
    dustrdz_neap=dustrdz[it_neap,:,:]
    dustrdz_spring=dustrdz[it_spring,:,:]

    #Make a figure for this section
    #get corners for pcolor plot
    X=xr.DataArray(np.concatenate((np.array([0]),ds['dd'].cumsum(dim='p').to_numpy())), dims='p')
    X=X-X.isel(p=-1)/2
    Y=np.concatenate((-h[np.newaxis,:], np.cumsum(dz,axis=0)-h[np.newaxis,:]),axis=0)
    # Y=(np.concatenate((Ydata[:,:,0,None],Ydata), axis=2) + np.concatenate((Ydata,Ydata[:,:,-1,None]), axis=2))/2
    Y=(np.concatenate((Y[:,np.newaxis,0],Y),axis=1) + np.concatenate((Y,Y[:,np.newaxis,-1]),axis=1))/2
    #Y=xr.DataArray(Ydata, coords={'time':ds.time}, dims=['time','z','p'])

    #plot
    vlim=np.max(np.abs(pfun.auto_lims(pg[it_neap,:,:])))
    fig, axs = plt.subplots(2, 3,figsize=(20,10),sharex=True)
    cs0=axs[0,0].pcolormesh(X,Y,pg_neap,cmap=cm.balance,norm=colors.CenteredNorm(halfrange=vlim)) #try twoslopenorm instead of colors.CenteredNorm()
    cs1=axs[0,1].pcolormesh(X,Y,dustrdz_neap,cmap=cm.balance,norm=colors.CenteredNorm(halfrange=vlim)) #try twoslopenorm instead of colors.CenteredNorm()
    cs2=axs[0,2].pcolormesh(X,Y,dustrdz_neap-pg_neap,cmap=cm.balance,norm=colors.CenteredNorm(halfrange=vlim)) #try twoslopenorm instead of colors.CenteredNorm()
    cs3=axs[1,0].pcolormesh(X,Y,pg_spring,cmap=cm.balance,norm=colors.CenteredNorm(halfrange=vlim)) #try twoslopenorm instead of colors.CenteredNorm()
    cs4=axs[1,1].pcolormesh(X,Y,dustrdz_spring,cmap=cm.balance,norm=colors.CenteredNorm(halfrange=vlim)) #try twoslopenorm instead of colors.CenteredNorm()
    cs5=axs[1,2].pcolormesh(X,Y,dustrdz_spring-pg_spring,cmap=cm.balance,norm=colors.CenteredNorm(halfrange=vlim)) #try twoslopenorm instead of colors.CenteredNorm()

    # cs0=axs[0,0].pcolormesh(X,Y,pg[it_neap,:,:],cmap=cm.balance,vmin=-1.75e-5,vmax=1.75e-5) #try twoslopenorm instead of colors.CenteredNorm()
    # cs1=axs[0,1].pcolormesh(X,Y,dustrdz[it_neap,:,:],cmap=cm.balance,vmin=-1.75e-5,vmax=1.75e-5) #try twoslopenorm instead of colors.CenteredNorm()
    # cs2=axs[0,2].pcolormesh(X,Y,pg[it_neap,:,:]-dustrdz[it_neap,:,:],cmap=cm.balance,vmin=-1.75e-5,vmax=1.75e-5) #try twoslopenorm instead of colors.CenteredNorm()
    # cs3=axs[1,0].pcolormesh(X,Y,pg[it_spring,:,:],cmap=cm.balance,vmin=-1.75e-5,vmax=1.75e-5) #try twoslopenorm instead of colors.CenteredNorm()
    # cs4=axs[1,1].pcolormesh(X,Y,dustrdz[it_spring,:,:],cmap=cm.balance,vmin=-1.75e-5,vmax=1.75e-5) #try twoslopenorm instead of colors.CenteredNorm()
    # cs5=axs[1,2].pcolormesh(X,Y,pg[it_spring,:,:]-dustrdz[it_spring,:,:],cmap=cm.balance,vmin=-1.75e-5,vmax=1.75e-5) #try twoslopenorm instead of colors.CenteredNorm()

    axs[0,0].set_title('PG neap')
    axs[0,1].set_title('Stress divergence neap')
    axs[0,2].set_title('Difference neap')
    axs[1,0].set_title('PG spring')
    axs[1,1].set_title('Stress divergence spring')
    axs[1,2].set_title('Difference spring')

    cb0=fig.colorbar(cs0, ax=axs[0,0])
    cb1=fig.colorbar(cs1, ax=axs[0,1])
    cb2=fig.colorbar(cs2, ax=axs[0,2])
    cb3=fig.colorbar(cs3, ax=axs[1,0])
    cb4=fig.colorbar(cs4, ax=axs[1,1])
    cb5=fig.colorbar(cs5, ax=axs[1,2])
    plt.suptitle(Ldir['gridname']+' '+sn)

    fig.tight_layout()
    #fig.savefig(out_dir / 'dsdx_spring_neap.png')
    fig.savefig(out_dir / 'pg_stressdiv_'+sn+'.png',dpi=300)

    #save fields for making multi-section plot
        # stz_dict[sn] = s_vs_z
    # dustrdztz_dict[sn] = dustrdz_vs_z
    zeta_dict[sn] = zeta_avg
    pg_center_neap_dict[sn]= pgcorr[it_neap,:,6]
    pg_center_spring_dict[sn]= pgcorr[it_spring,:,6]
    dustrdz_center_neap_dict[sn]= dustrdz[it_neap,:,6]
    dustrdz_center_spring_dict[sn]= dustrdz[it_spring,:,6]
    zr_center_dict[sn] = zr[:,6]
#fig = plt.figure(figsize=(20,15))
#gs = fig.add_gridspec(nrows=2,ncols=3,width_ratios=[1,1,1],height_ratios=[1,1])

# Make a plot with all the section at the central water column
fig, [ax0,ax1] = plt.subplots(1, 2,figsize=(20,10),sharex=True)
ci=0
for sn in sect_list:
    zrc=zr_center_dict[sn]
    pgcneap=pg_center_neap_dict[sn]
    pgcspring=pg_center_spring_dict[sn]
    dustrdzcneap=dustrdz_center_neap_dict[sn]
    dustrdzcspring=dustrdz_center_spring_dict[sn]
    ax0.plot(-pgcneap,zrc,color=color_list1[ci],label=sn+' PG')
    ax0.plot(dustrdzcneap,zrc,color=color_list2[ci],label=sn+' stress div')
    ax0.plot(dustrdzcneap-pgcneap,zrc,ls=':',color=color_list1[ci],label=sn+' difference')
    ax1.plot(-pgcspring,zrc,color=color_list1[ci],label=sn+' PG')
    ax1.plot(dustrdzcspring,zrc,color=color_list2[ci],label=sn+' stress div')
    ax1.plot(dustrdzcspring-pgcspring,zrc,color=color_list1[ci],label=sn+' difference')
    ax0.set_title('Neap')
    ax1.set_title('Spring')
    ci=ci+1

fig.tight_layout()
#fig.savefig(out_dir / 'dsdx_spring_neap.png')
fig.savefig(out_dir / 'pg_stressdiv_summary.png',dpi=300)

# # get dx for ds/dx #MIGHT CHANGE THIS FOR MORE PAIRS ALONG THE ESTUARY
# dxlist = []
# for i in range(len(sect_list)-1):
#     dx, ang = sw.dist([lat_dict[sect_list[i]],lat_dict[sect_list[i+1]]],
#         [lon_dict[sect_list[i]],lon_dict[sect_list[i+1]]],
#         units='km')
#     dxlist.append(dx[0])
    
# # z for plotting
# z = bin_edges[:-1] + np.diff(bin_edges)/2
# z2 = bin_edges2[:-1] + np.diff(bin_edges2)/2 #use with dustrdzf

# # trim to only use overlapping z range
# mask = z == z
# Stz_dict = dict() # trimmed version of stz
# # mask is initialized as all True
# for sn in sect_list:
#     s0z = stz_dict[sn][0,:]
#     mask = mask & ~np.isnan(s0z)
# for sn in sect_list:
#     stz = stz_dict[sn]
#     Stz_dict[sn] = stz[:,mask]
# Z = z[mask] # trimmed version of z

# Sz_dict = dict() # time-mean of each section s(z)
# St_dict = dict() # depth-mean of each section s(t)
# for sn in sect_list:
#     # Stz is a trimmed array of s(t,z), daily
#     Stz = Stz_dict[sn]
#     Sz_dict[sn] = np.mean(Stz,axis=0)
#     St_dict[sn] = np.nanmean(Stz,axis=1)

# # trim ustr to only use overlapping z range
# mask2 = z2 == z2
# DUSTRDZtz_dict = dict() # trimmed version of ustrtz
# # mask is initialized as all True
# for sn in sect_list:
#     dustrdz0z = dustrdztz_dict[sn][0,:]
#     mask2 = mask2 & ~np.isnan(dustrdz0z)
# for sn in sect_list:
#     dustrdztz = dustrdztz_dict[sn]
#     DUSTRDZtz_dict[sn] = dustrdztz[:,mask2]
# Z2 = z2[mask2] # trimmed version of z

# # vertical derivative of ustr #NOW WE ARE CALCULATING THIS BEFORE BINNING
# dustrdz_dict = dict() #vertical derivative of u stress
# for sn in sect_list:
#     # USTRtz is a trimmed array of s(t,z), daily
#     USTRtz = USTRtz_dict[sn]
#     dustrdz_dict[sn] = np.diff(USTRtz,axis=1)/np.diff(Z2)

# # useful time vectors
# dti = pd.DatetimeIndex(otdt)
# yd = dti.dayofyear
# year = otdt[0].year
    
# # plotting
# plt.close('all')
# pfun.start_plot(figsize=(12,12))
# fig = plt.figure()

# #c_list = ['m','r','orange','g','b','violet']
# c_list = ['tab:red','tab:orange','tab:green','tab:cyan','tab:blue'] #COLORS FOR SHORT SECTION LIST ON SILL
# c_dict = dict(zip(sect_list,c_list))

# # map
# ax0 = fig.add_subplot(321)
# lon0 = lon_vec_dict[sect_list[0]]
# lat0 = lat_vec_dict[sect_list[0]]
# lon1 = lon_vec_dict[sect_list[-1]]
# lat1 = lat_vec_dict[sect_list[-1]]
# lonmin = np.min(np.concatenate((lon0,lon1)))
# lonmax = np.max(np.concatenate((lon0,lon1)))
# latmin = np.min(np.concatenate((lat0,lat1)))
# latmax = np.max(np.concatenate((lat0,lat1)))
# # ax.pcolormesh() #need to add something here
# cs = ax0.pcolormesh(plon, plat, zm, vmin=-300, vmax=20, cmap='Greys_r')
# for sn in sect_list:
#     ax0.plot(lon_vec_dict[sn], lat_vec_dict[sn], '.',color=c_dict[sn])
# #pfun.add_coast(ax,color='gray',linewidth=2) #add coast doesn't work for idealized model
# mpad = .2
# # ax0.axis([lonmin-mpad, lonmax+mpad, latmin-mpad, latmax+mpad])
# ax0.axis([lonmin-mpad, lonmax+mpad, 44.9, 45.1])
# pfun.dar(ax0)
# ax0.text(.05,.9,'(a) Section Locations',color='k',fontweight='bold',
#     transform=ax0.transAxes,bbox=pfun.bbox)
# ax0.set_xlabel('Longitude')
# ax0.set_ylabel('Latitude')

# ax1 = fig.add_subplot(322)
# if False:
#     for sn in sect_list:
#         ax.plot(Sz_dict[sn],Z,'-',color=c_dict[sn])
#     ax.text(.05,.1,'(b) Time-Mean S(z)',color='k',fontweight='bold',
#         transform=ax.transAxes,bbox=pfun.bbox)
# else:
#     # selected spring and neap; hard coded for 2018 Admiralty Inlet #NEED TO CHANGE THIS
#     # it_neap = zfun.find_nearest_ind(yd,233)
#     # it_spring = zfun.find_nearest_ind(yd,253)
#     # it_neap = zfun.find_nearest_ind(yd,279)
#     # it_spring = zfun.find_nearest_ind(yd,287)
#     it_neap = zfun.find_nearest_ind(dti,TE['t_neap'])
#     it_spring = zfun.find_nearest_ind(dti,TE['t_spring'])
#     for sn in sect_list:
#         ax1.plot(Stz_dict[sn][it_neap,:],Z,'-',color=c_dict[sn])
#         ax1.plot(Stz_dict[sn][it_spring,:],Z,'--',color=c_dict[sn])
#     ax1.text(.05,.1,'(b) Neap and Spring S(z)',color='k',fontweight='bold',
#         transform=ax1.transAxes,bbox=pfun.bbox)
# ax1.set_xlabel('Salinity')
# ax1.set_ylabel('Z [m]')
# ax1.grid(True)


# it_neap = zfun.find_nearest_ind(dti,TE['t_neap'])
# it_spring = zfun.find_nearest_ind(dti,TE['t_spring'])

# ax1 = fig.add_subplot(121)
# for sn in sect_list:
# #    ax1.plot(DUSTRDZtz_dict[sn][it_neap,:],Z2,'-',color=c_dict[sn])
#     ax1.plot(dudz_hourly[1,:,6],zw[1:-1,6],'-',color=c_dict[sn])
# ax1.text(.05,.1,'Neap u stress',color='k',fontweight='bold',transform=ax1.transAxes,bbox=pfun.bbox)
# ax1.set_xlabel('dudz center b5')
# ax1.set_ylabel('Z [m]')
# ax1.grid(True)

# ax2 = fig.add_subplot(122)
# for sn in sect_list:
#     ax2.plot(u_hourly[1,:,6],zr[:,6],'--',color=c_dict[sn])
# ax2.text(.05,.1,'Spring u stress',color='k',fontweight='bold',transform=ax2.transAxes,bbox=pfun.bbox)
# ax2.set_xlabel('u center b5')
# ax2.set_ylabel('Z [m]')
# ax2.grid(True)


# ax2 = fig.add_subplot(312)
# for sn in sect_list:
#     ax2.plot(dti,St_dict[sn],'-',color=c_dict[sn])
# ax2.text(.05,.9,'(c) Depth-Mean S(t)',color='k',fontweight='bold',
#     transform=ax2.transAxes,bbox=pfun.bbox)
# #ax.set_xlim(0,365)
# #ax2.set_xlim(246,365) #change this to monthday or something!!
# # ax.set_xlabel('Yearday ' + str(year))
# #ax2.grid(axis='x')
# if True:
#     ax2.axvline(x=dti[it_neap],linestyle='-',color='gray',linewidth=2)
#     ax2.axvline(x=dti[it_spring],linestyle='--',color='gray',linewidth=2)

# #ax3a = fig.add_subplot(313)
# ax3a = fig.add_subplot(111)
# dti = pd.DatetimeIndex(otdt)
# yd = dti.dayofyear
# year = otdt[0].year
# ##ax.plot(yd,St_dict[sect_list[0]]-St_dict[sect_list[-1]],'-',color='k')
# # ax.plot(yd,St_dict[sect_list[0]]-St_dict[sect_list[2]],'-',color=c_dict[sect_list[0]])
# # ax.plot(yd,St_dict[sect_list[-3]]-St_dict[sect_list[-1]],'-',color=c_dict[sect_list[-1]])
# ##ax.text(.05,.9,'(d) Total Along-Section Change in Depth-Mean Salinity',color='k',fontweight='bold',transform=ax.transAxes,bbox=pfun.bbox)
# for i in range(len(sect_list)-1):
#     ax3a.plot(dti,(St_dict[sect_list[i]]-St_dict[sect_list[i+1]])/dxlist[i],'-',color=c_dict[sect_list[i]]) #PLOT ds/dx INSTEAD OF SALINITY CHANGE
#     ax3a.plot(dti,(St_dict[sect_list[i]]-St_dict[sect_list[i+1]])/dxlist[i],':',color=c_dict[sect_list[i+1]]) #PLOT AGAIN IN SECOND COLOR TO MAKE TWO COLOR DASHED LINE
# ax3a.text(.05,.05,r'(d) $\partial S/\partial x\ [g\ kg^{-1}\ km^{-1}]$ between pairs of sections',color='k',fontweight='bold',transform=ax3a.transAxes,bbox=pfun.bbox)
# #ax.set_xlim(0,365)
# ax3a.set_xlim(246,365) #change this to monthday or something!!
# ax3a.set_xlabel('Yearday ' + str(year))
# #ax3a.grid(axis='x')
# # add Qprism
# # ax2b = ax2.twinx()
# ax3b = ax3a.twinx()
# Qprism_sectavg = 0.5*(Qprism_dict[sect_list[0]]+Qprism_dict[sect_list[-1]])/1000
# # ax2b.plot(dti,Qprism_sectavg,'-',color='tab:purple',linewidth=3,alpha=.4)
# # ax2b.set_ylim(bottom=0)
# ax3b.plot(dti,Qprism_sectavg,'-',color='tab:purple',linewidth=3,alpha=.4)
# ax3b.set_ylim(bottom=0)

# snmid=(np.max(Qprism_sectavg)+np.min(Qprism_sectavg))/2
# snbg=np.where(Qprism_sectavg>snmid, 1, 0)
# ax3b.pcolor(dti, ax3b.get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True) #add grey bars for qprism
# # ax2.pcolor(dti, ax2.get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True) #also add to subplot above grey bars for qprism
# ax3a.set_xlim(pd.Timestamp('2020-10-01'), pd.Timestamp('2020-10-31'))
# ax3a.xaxis.set_major_formatter(mdates.DateFormatter('%-d'))
# ax3a.set_xlabel('Day')
# # ax2.set_xlim(pd.Timestamp('2020-10-01'), pd.Timestamp('2020-10-31'))
# # ax2.xaxis.set_major_formatter(mdates.DateFormatter('%-d'))
# # ax2.set_xlabel('Day')
# # ax2.grid(False)
# # ax3a.grid(False)
# # ax2.grid(axis='y')
# ax3a.grid(axis='y')
# if True:
#     ax3a.axvline(x=dti[it_neap],linestyle='-',color='gray',linewidth=2)
#     ax3a.axvline(x=dti[it_spring],linestyle='--',color='gray',linewidth=2)

# # ax2b.text(.95,.9,r'$Q_{prism}\ [10^{3}m^{3}s^{-1}]$', color='tab:purple', 
# # transform=ax2.transAxes, ha='right',
# # bbox=pfun.bbox)
# # ax2b.xaxis.label.set_color('tab:purple')
# # ax2b.tick_params(axis='y', colors='tab:purple')

# ax3b.text(.95,.9,r'$Q_{prism}\ [10^{3}m^{3}s^{-1}]$', color='tab:purple', 
#     transform=ax3a.transAxes, ha='right',
#     bbox=pfun.bbox)
# ax3b.xaxis.label.set_color('tab:purple')
# ax3b.tick_params(axis='y', colors='tab:purple')
# #ax.set_xlim(0,365)
# #ax3a.set_xlim(246,365) #change this to monthday or something!!
# # ax3a.set_ylim(bottom=0)
# if True:
#     ax3a.axvline(x=yd[it_neap],linestyle='-',color='gray',linewidth=2)
#     ax3a.axvline(x=yd[it_spring],linestyle='--',color='gray',linewidth=2)
    
# fig.tight_layout()
# #fig.savefig(out_dir / 'dsdx_spring_neap.png')
# fig.savefig(out_dir / 'pg_stressdiv.png',dpi=300)

# plt.show()