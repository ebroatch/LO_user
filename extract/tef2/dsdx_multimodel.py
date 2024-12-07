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
import tef_fun

from lo_tools import Lfun, zrfun, zfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

# output location
out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
out_dir = out_dir0 / ('dsdx_plots_' + Ldir['ds0'] + '_' + Ldir['ds1'])
#out_dir = Ldir['parent'] / 'LPM_output' / 'extract'/ 'tef_exdyn'
Lfun.make_dir(out_dir)

# grids to use
# gridnamelist = ['sill5km','sill10km','sill20kmdeep','sill40km','sill80km']
# taglist = ['t0','t2','t2','t2','t2']
# gtaglist = ['sill5km_t0','sill10km_t2','sill20kmdeep_t2','sill40km_t2','sill80km_t2']
# gtagexlist = ['sill5km_t0_xa0','sill10km_t2_xa0','sill20kmdeep_t2_xa0','sill40km_t2_xa0','sill80km_t2_xa0']
# gridlist = [Ldir['data'] / 'grids' / 'sill5km', Ldir['data'] / 'grids' / 'sill10km', Ldir['data'] / 'grids' / 'sill20kmdeep', Ldir['data'] / 'grids' / 'sill40km', Ldir['data'] / 'grids' / 'sill80km']
# silllenstr = ['5km','10km','20km','40km','80km']
# silllen = [5e3,10e3,20e3,40e3,80e3]
gridnamelist = ['sill5km','sill40km']
taglist = ['t0','t2']
gtaglist = ['sill5km_t0','sill40km_t2']
gtagexlist = ['sill5km_t0_xa0','sill40km_t2_xa0']
gridlist = [Ldir['data'] / 'grids' / 'sill5km', Ldir['data'] / 'grids' / 'sill40km']
silllenstr = ['5km','40km']
silllen = [5e3,40e3]

# Define sections to work on.
# Generally choose [seaward, landward]
sect_list = ['b1','b2','b3','b4','b5'] #SHORT LIST SILL ONLY
#sect_list = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5'] #FULL LIST
# plot_color = ['tab:cyan',plt.cm.Dark2(0),'tab:olive',plt.cm.Dark2(5),'tab:brown']
plot_color = ['tab:cyan',plt.cm.Dark2(0),plt.cm.Dark2(5),'tab:brown'] #need one less color for the pairs

# plotting
plt.close('all')
pfun.start_plot(figsize=(8,8),fs=14)
fig, axs = plt.subplots(3, 1, sharex=True,figsize=(8,8),gridspec_kw={'height_ratios': [1,4,4]}) #change # rows if adding more models

for gi in range(len(gtagexlist)):
    Ldir['gridname']=gridnamelist[gi]
    Ldir['tag']=taglist[gi]
    Ldir['gtag']=gtaglist[gi]
    Ldir['gtagex']=gtagexlist[gi]
    Ldir['grid']=gridlist[gi]

    # gctag and location of tef2 section definitions
    gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
    tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

    # get sect_df with the section point locations
    sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
    sect_df = pd.read_pickle(sect_df_fn)

    # where to find the extracted sections
    in_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
    in_dir = in_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
    in_dir2 = in_dir0 / ('bulk_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])

    # get tide info from the tide excursion calculator
    excur_dir = in_dir0 / ('tide_excursion_' + Ldir['ds0'] + '_' + Ldir['ds1'])
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

    # make vn_list by inspecting the first section
    ds = xr.open_dataset(in_dir / (sect_list[0] + '.nc'))
    vn_list = [item for item in ds.data_vars \
        if (len(ds[item].dims) == 3) and (item not in ['vel','DZ'])]
    ds.close()

    print('\nCalculating s(z) on mutiple sections')
    print(str(in_dir))

    tt00 = time()

    stz_dict = dict()
    lon_dict = dict()
    lat_dict = dict()
    lon_vec_dict = dict()
    lat_vec_dict = dict()
    H_dict = dict()
    Qprism_dict = dict()

    for sn in sect_list:
        tt0 = time()
        print(sn)

        # load fields
        ds = xr.open_dataset(in_dir / (sn + '.nc'))
        salt_hourly = ds.salt.values
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
        # what if there is only one p?
        H_dict[sn] = h.mean()
        
        # get section area (ssh=0)
        dz = np.diff(zw,axis=0)
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
        
        # Then we want to form a time series of s(z)
        NZ = 100
        s_vs_z = np.nan * np.ones((NT,NZ))
        for tt in range(NT):
            sf = salt[tt,:,:].squeeze().flatten()
            # scipy.stats.binned_statistic(x, values, statistic='mean', bins=10, range=None)
            ##bs = binned_statistic(zf, sf, statistic='mean', bins=NZ, range=(-500,0))
            bs = binned_statistic(zf, sf, statistic='mean', bins=NZ, range=(-200,0)) #change to -200 max depth in estuary
            s_vs_z[tt,:] = bs.statistic
        bin_edges = bs.bin_edges
        
        ot = ds['time'].to_numpy()
        # do a little massaging of ot
        dti = pd.to_datetime(ot) # a pandas DatetimeIndex with dtype='datetime64[ns]'
        dt = dti.to_pydatetime() # an array of datetimes
        ot = np.array([Lfun.datetime_to_modtime(item) for item in dt])
        ot = ot[pad:-pad+1:24]
        # also make an array of datetimes to save as the ot variable
        otdt = np.array([Lfun.modtime_to_datetime(item) for item in ot])
        
        stz_dict[sn] = s_vs_z

    # get dx for ds/dx #MIGHT CHANGE THIS FOR MORE PAIRS ALONG THE ESTUARY
    dxlist = []
    for i in range(len(sect_list)-1):
        dx, ang = sw.dist([lat_dict[sect_list[i]],lat_dict[sect_list[i+1]]],
            [lon_dict[sect_list[i]],lon_dict[sect_list[i+1]]],
            units='km')
        dxlist.append(dx[0])
        
    # z for plotting
    z = bin_edges[:-1] + np.diff(bin_edges)/2

    # trim to only use overlapping z range
    mask = z == z
    Stz_dict = dict() # trimmed version of stz
    # mask is initialized as all True
    for sn in sect_list:
        s0z = stz_dict[sn][0,:]
        mask = mask & ~np.isnan(s0z)
    for sn in sect_list:
        stz = stz_dict[sn]
        Stz_dict[sn] = stz[:,mask]
    Z = z[mask] # trimmed version of z
        
    Sz_dict = dict() # time-mean of each section s(z)
    St_dict = dict() # depth-mean of each section s(t)
    for sn in sect_list:
        # Stz is a trimmed array of s(t,z), daily
        Stz = Stz_dict[sn]
        Sz_dict[sn] = np.mean(Stz,axis=0)
        St_dict[sn] = np.nanmean(Stz,axis=1)

    # useful time vectors
    dti = pd.DatetimeIndex(otdt)
    yd = dti.dayofyear
    year = otdt[0].year

    ##ax.text(.05,.9,'(d) Total Along-Section Change in Depth-Mean Salinity',color='k',fontweight='bold',transform=ax.transAxes,bbox=pfun.bbox)
    # ax3a.text(.05,.05,r'(d) $\partial S/\partial x\ [g\ kg^{-1}\ km^{-1}]$ between pairs of sections',color='k',fontweight='bold',transform=ax3a.transAxes,bbox=pfun.bbox)

    for i in range(len(sect_list)-1):
        axs[gi+1].plot(dti,(St_dict[sect_list[i+1]]-St_dict[sect_list[i]])/dxlist[i],'-',color=plot_color[i],label=sect_list[i+1]+'-'+sect_list[i]) #fix sign (match with label and coord system)
        # axs[gi+1].plot(dti,(St_dict[sect_list[i]]-St_dict[sect_list[i+1]])/dxlist[i],'-',color=plot_color[i],label=sect_list[i+1]+'-'+sect_list[i]) #PLOT ds/dx INSTEAD OF SALINITY CHANGE
        # axs[gi+1].plot(dti,(St_dict[sect_list[i]]-St_dict[sect_list[i+1]])/dxlist[i],':',color=c_dict[sect_list[i+1]]) #PLOT AGAIN IN SECOND COLOR TO MAKE TWO COLOR DASHED LINE
    axs[gi+1].plot(dti,(St_dict[sect_list[4]]-St_dict[sect_list[0]])/(dxlist[0]+dxlist[1]+dxlist[2]+dxlist[3]),'--',color='k',label=sect_list[4]+'-'+sect_list[0]) #plot ds/dx across whole sill
    
    if gi==0:
        #Add Qprism and grey bars
        axs[0].set_ylim(20,80)
        axs[0].set_yticks([20,50,80])
        # axs[1].set_ylim(-0.1,0.5) #to match when placed side by side
        # axs[2].set_ylim(0.02,0.14)
        axs[1].set_ylim(-0.5,0.1) #to match when placed side by side
        axs[1].axhline(y=0,color='tab:gray',linewidth=2)
        axs[2].set_ylim(-0.14,-0.02)

        sect_name='b3'
        bulk = xr.open_dataset(in_dir2 / (sect_name + '.nc'))
        tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir2, sect_name)
        tef_df['Q_p'] = tef_df['q_p']/1000
        tef_df['Q_m'] = tef_df['q_m']/1000
        tef_df['Q_prism']=tef_df['qprism']/1000
        ot = bulk.time.values
        lw=2
        axs[0].plot(ot,tef_df['Q_prism'].to_numpy(), color='tab:gray', linewidth=lw)
        snmid=(np.max(tef_df['Q_prism'].loc['2020-10-01':'2020-10-31'])+np.min(tef_df['Q_prism'].loc['2020-10-01':'2020-10-31']))/2
        snbg=np.where(tef_df['Q_prism'].to_numpy()>snmid, 1, 0)
        axs[0].pcolor(ot, axs[0].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
        axs[1].pcolor(ot, axs[1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
        axs[2].pcolor(ot, axs[2].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
        

axs[2].set_xlim(pd.Timestamp('2020-10-01'), pd.Timestamp('2020-10-31'))
axs[2].xaxis.set_major_formatter(mdates.DateFormatter('%j'))
axs[2].set_xlabel('Yearday')
axs[1].legend(loc='lower right',fontsize=12)

axs[0].grid(True)
axs[1].grid(True)
axs[2].grid(True)

axs[0].text(.02, .8, 'A', horizontalalignment='left', verticalalignment='top', transform=axs[0].transAxes, fontsize=14, fontweight='bold')
axs[1].text(.02, .95, 'B', horizontalalignment='left', verticalalignment='top', transform=axs[1].transAxes, fontsize=14, fontweight='bold')
axs[2].text(.02, .95, 'C', horizontalalignment='left', verticalalignment='top', transform=axs[2].transAxes, fontsize=14, fontweight='bold')
axs[1].text(0.99,0.98,'5km sill model',fontsize=10,ha='right',va='top',transform=axs[1].transAxes)
axs[2].text(0.99,0.98,'40km sill model',fontsize=10,ha='right',va='top',transform=axs[2].transAxes)

axs[0].set_ylabel('$Q_{prism}$\n(5km b3)\n$[10^{3}\ m^{3}s^{-1}]$')
axs[1].set_ylabel('$ds/dx$\n$[g\ kg^{-1}km^{-1}]$')
axs[2].set_ylabel('$ds/dx$\n$[g\ kg^{-1}km^{-1}]$')

fig.tight_layout()
#fig.savefig(out_dir / 'dsdx_spring_neap.png')
fig.savefig(out_dir / 'dsdx_spring_neap_pairs_multimodel.png',dpi=300)

plt.show()