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
from cmocean import cm

from lo_tools import Lfun, zrfun, zfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing
# use the arg -sect_name to pick section for plots

# import tef_fun


gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
in_dir2 = out_dir0 / ('processed_' + Ldir['ds0'] + '_' + Ldir['ds1']) #in_dir for processed to get qnet
out_dir = out_dir0 / ('spatial_structure_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

# sect_list = [item.name for item in in_dir.glob('*.nc')]
# if Ldir['testing']:
#     sect_list = ['jdf3.nc']
sect_list = [Ldir['sect_name']+'.nc'] #section to plot #not sure if this syntax is correct
    
# make vn_list by inspecting the first section
ds = xr.open_dataset(in_dir / sect_list[0])
vn_list = [item for item in ds.data_vars \
    if (len(ds[item].dims) == 3) and (item not in ['vel','DZ'])]
ds.close()

print('\nProcessing standard decomposition:')
#print(str(in_dir))

tt00 = time()
pad=36

t_spring_ebb = pd.Timestamp('2020-07-01 04:00:00')
t_spring_flood = pd.Timestamp('2020-07-01 10:00:00')
t_neap_ebb = pd.Timestamp('2020-07-08 10:00:00')
t_neap_flood = pd.Timestamp('2020-07-08 16:00:00')


for ext_fn in sect_list:
    tt0 = time()
    print(ext_fn)

    # name output file
    out_fn = ext_fn.replace('.nc','.p')

    # load fields
    ds = xr.open_dataset(in_dir / ext_fn)
    # V = dict()
    # for vn in vn_list:
    #     V[vn] = ds[vn].to_numpy()
    # V['salt2'] = V['salt']*V['salt']

    ot = ds['time'].to_numpy() #shape NT
    zeta = ds['zeta'].to_numpy() #shape NT,NX
    h = ds['h'].to_numpy() #shape NX
    dd = ds['dd'].to_numpy() #shape NX
    dz = ds['DZ'].to_numpy() #shape NT,NZ,NX
    u = ds['vel'].to_numpy() #shape NT,NZ,NX
    s = ds['salt'].to_numpy() #shape NT,NZ,NX
    
    dA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() #shape NT,NZ,NX
    H = np.sum(ds['DZ'].to_numpy(),axis=1) #shape NT,NX
    q = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['vel'].to_numpy() #shape NT,NZ,NX
    sdA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['salt'].to_numpy() #shape NT,NZ,NX
    # V['q'] = q

    #fig, axs = plt.subplots(2, 2,figsize=(15,15))
    fig = plt.figure(figsize=(20,15))
    gs = fig.add_gridspec(nrows=3,ncols=4,width_ratios=[1,1,1,1],height_ratios=[2,2,1])
    #gs = fig.add_gridspec(nrows=5,ncols=4,width_ratios=[1,1,1,1],height_ratios=[2,2,1,1,1])
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[0,2])
    ax4 = fig.add_subplot(gs[0,3])
    ax5 = fig.add_subplot(gs[1,0])
    ax6 = fig.add_subplot(gs[1,1])
    ax7 = fig.add_subplot(gs[1,2])
    ax8 = fig.add_subplot(gs[1,3])
    ax9 = fig.add_subplot(gs[2,:])
    # ax10 = fig.add_subplot(gs[3,:])
    # ax11 = fig.add_subplot(gs[4,:])

    X=xr.DataArray(np.concatenate((np.array([0]),ds['dd'].cumsum(dim='p').to_numpy())), dims='p')
    X=X-X.isel(p=-1)/2
    Ydata=np.concatenate((np.expand_dims((xr.ones_like(ds['DZ'].isel(z=0))*(-ds['h'])).to_numpy(), axis=1) , (ds['DZ'].cumsum(dim='z')-ds['h']).to_numpy()),axis=1)
    Ydata=(np.concatenate((Ydata[:,:,0,None],Ydata), axis=2) + np.concatenate((Ydata,Ydata[:,:,-1,None]), axis=2))/2
    Y=xr.DataArray(Ydata, coords={'time':ds.time}, dims=['time','z','p'])
    
    ulim=0.8
    slimmin=20
    slimmax=34
    cs1=ax1.pcolormesh(X,Y.sel(time=t_spring_ebb),ds['vel'].sel(time=t_spring_ebb),cmap=cm.balance,vmin=-ulim,vmax=ulim)
    cs2=ax2.pcolormesh(X,Y.sel(time=t_spring_flood),ds['vel'].sel(time=t_spring_flood),cmap=cm.balance,vmin=-ulim,vmax=ulim)
    cs3=ax3.pcolormesh(X,Y.sel(time=t_neap_ebb),ds['vel'].sel(time=t_neap_ebb),cmap=cm.balance,vmin=-ulim,vmax=ulim)
    cs4=ax4.pcolormesh(X,Y.sel(time=t_neap_flood),ds['vel'].sel(time=t_neap_flood),cmap=cm.balance,vmin=-ulim,vmax=ulim)
    cs5=ax5.pcolormesh(X,Y.sel(time=t_spring_ebb),ds['salt'].sel(time=t_spring_ebb),cmap=cm.haline,vmin=slimmin,vmax=slimmax)
    cs6=ax6.pcolormesh(X,Y.sel(time=t_spring_flood),ds['salt'].sel(time=t_spring_flood),cmap=cm.haline,vmin=slimmin,vmax=slimmax)
    cs7=ax7.pcolormesh(X,Y.sel(time=t_neap_ebb),ds['salt'].sel(time=t_neap_ebb),cmap=cm.haline,vmin=slimmin,vmax=slimmax)
    cs8=ax8.pcolormesh(X,Y.sel(time=t_neap_flood),ds['salt'].sel(time=t_neap_flood),cmap=cm.haline,vmin=slimmin,vmax=slimmax)

    fig.colorbar(cs1, ax=ax1)
    fig.colorbar(cs2, ax=ax2)
    fig.colorbar(cs3, ax=ax3)
    fig.colorbar(cs4, ax=ax4)
    fig.colorbar(cs5, ax=ax5)
    fig.colorbar(cs6, ax=ax6)
    fig.colorbar(cs7, ax=ax7)
    fig.colorbar(cs8, ax=ax8)

    ax1.set_title('u spring ebb', c='tab:olive', fontweight='bold')
    ax2.set_title('u spring flood', c='tab:green', fontweight='bold')
    ax3.set_title('u neap ebb', c='tab:purple', fontweight='bold')
    ax4.set_title('u neap flood', c='tab:blue', fontweight='bold')
    ax5.set_title('salt spring ebb', c='tab:olive', fontweight='bold')
    ax6.set_title('salt spring flood', c='tab:green', fontweight='bold')
    ax7.set_title('salt neap ebb', c='tab:purple', fontweight='bold')
    ax8.set_title('salt neap flood', c='tab:blue', fontweight='bold')

    # load fields
    ds2 = xr.open_dataset(in_dir2 / ext_fn)
    ax9.plot(ds2['time'],ds2['qnet'])
    ax9.axvline(x=t_spring_ebb, c='tab:olive', linewidth=3)
    ax9.axvline(x=t_spring_flood, c='tab:green', linewidth=3)
    ax9.axvline(x=t_neap_ebb, c='tab:purple', linewidth=3)
    ax9.axvline(x=t_neap_flood, c='tab:blue', linewidth=3)  
    ax9.grid(True)
    ax9.set_xlim(pd.Timestamp('2020-06-25'), pd.Timestamp('2020-07-10')) #to see tidal cycle zoom
    # ax9.set_xlim(pd.Timestamp('2020-06-30'), pd.Timestamp('2020-07-02')) #to see tidal cycle zoom
    ax9.set_title('qnet tidal transport')

    # ax10.plot(ds2['time'],ds2['qnet'])
    # ax10.axvline(x=t_spring_flood, c='tab:green')
    # ax10.axvline(x=t_spring_ebb, c='tab:olive')
    # ax10.axvline(x=t_neap, c='tab:blue')
    # ax10.axvline(x=t_neap, c='tab:purple')
    # ax10.grid(True)
    # ax10.set_xlim(pd.Timestamp('2020-07-08'), pd.Timestamp('2020-07-09')) #to see tidal cycle zoom
    # ax10.set_title('qnet tidal transport')

    # qprism=zfun.lowpass(np.abs(ds2['qnet'].values-zfun.lowpass(ds2['qnet'].values, f='godin',nanpad=False)), f='godin')/2
    # ax11.plot(ds2['time'],qprism)
    # ax11.axvline(x=t_spring_flood, c='tab:green')
    # ax11.axvline(x=t_spring_ebb, c='tab:olive')
    # ax11.axvline(x=t_neap, c='tab:blue')
    # ax11.axvline(x=t_neap, c='tab:purple')
    # ax11.grid(True)
    # ax11.set_xlim(pd.Timestamp('2020-07-08T12'), pd.Timestamp('2020-07-08T16')) #to see tidal cycle zoom
    # ax11.set_ylim(10850,10950)
    # ax11.set_title('qprism')
    
    fig.suptitle(Ldir['sect_name'])
    plt.savefig(out_dir / (Ldir['sect_name'] + '.png'))
    ds.close()
    ds2.close()
    NT, NZ, NX = q.shape



    # A = np.sum(dA, axis=(1,2)) #shape NT
    # A0 = zfun.lowpass(A, f='godin')[pad:-pad+1] #shape NT-72
    # dA0 = zfun.lowpass(dA, f='godin')[pad:-pad+1, :, :] #shape NT-72,NZ,NX

    # u0 = (zfun.lowpass(np.sum(q, axis=(1,2)), f='godin')[pad:-pad+1])/A0 #shape NT-72
    # s0 = (zfun.lowpass(np.sum(sdA, axis=(1,2)), f='godin')[pad:-pad+1])/A0 #shape NT-72

    # u1 = (zfun.lowpass(q, f='godin')[pad:-pad+1, :, :])/dA0 - u0[:, np.newaxis, np.newaxis] #shape NT-72,NZ,NX newaxis for broadcasting
    # s1 = (zfun.lowpass(sdA, f='godin')[pad:-pad+1, :, :])/dA0 - s0[:, np.newaxis, np.newaxis] #shape NT-72,NZ,NX newaxis for broadcasting

    # u2 = u[pad:-pad+1, :, :] - u1 - u0[:, np.newaxis, np.newaxis] #shape NT-72,NZ,NX
    # s2 = s[pad:-pad+1, :, :] - s1 - s0[:, np.newaxis, np.newaxis] #shape NT-72,NZ,NX
    # dA2 = dA[pad:-pad+1, :, :] #shape NT-72,NZ,NX

    # u2L = (np.sum(u2*dz[pad:-pad+1, :, :],axis=1))/(H[pad:-pad+1, :]) #shape NT-72,NX
    # u2V = u2 - np.expand_dims(u2L,axis=1) #shape NT-72,NZ,NX expand_dims for broadcasting
    # s2L = (np.sum(s2*dz[pad:-pad+1, :, :],axis=1))/(H[pad:-pad+1, :]) #shape NT-72,NX
    # s2V = s2 - np.expand_dims(s2L,axis=1) #shape NT-72,NZ,NX expand_dims for broadcasting

    # ssh = np.mean(zeta, axis=1)[pad:-pad+1] #shape NT-72
    # ssh_lp = zfun.lowpass(np.mean(zeta, axis=1), f='godin')[pad:-pad+1] #shape NT-72
    
    # FR = (u0*s0*A0) #shape NT-72
    # FE = np.sum(u1*s1*dA0, axis=(1,2)) #shape NT-72
    # FT = zfun.lowpass(np.sum(u2*s2*dA2, axis=(1,2)), f='godin')[pad:-pad+1] #shape NT-144

    # FTL = zfun.lowpass(np.sum(u2L*s2L*H[pad:-pad+1, :]*dd, axis=1), f='godin')[pad:-pad+1] #shape NT-144
    # FTV = zfun.lowpass(np.sum(u2V*s2V*dA[pad:-pad+1, :, :], axis=(1,2)), f='godin')[pad:-pad+1] #shape NT-144

    # F = zfun.lowpass(np.sum(u*s*dA, axis=(1,2)), f='godin')[pad:-pad+1] #shape NT-72

    # SD = dict()
    # SD['u0']=u0[pad:-pad+1]
    # SD['s0']=s0[pad:-pad+1]
    # SD['A0']=A0[pad:-pad+1]
    # SD['FR']=FR[pad:-pad+1]
    # SD['FE']=FE[pad:-pad+1]
    # SD['FT']=FT
    # SD['FTL']=FTL
    # SD['FTV']=FTV
    # SD['F']=F[pad:-pad+1]
    # SD['ssh']=ssh[pad:-pad+1] #fix typo
    # SD['ssh_lp']=ssh_lp[pad:-pad+1]
    # SD['ot']=(ot[pad:-pad+1])[pad:-pad+1]
    # #pickle.dump(SD, open(out_dir / out_fn, 'wb'))
    print('  elapsed time for section = %d seconds' % (time()-tt0))
    sys.stdout.flush()


#print('\nTotal elapsed time for standard decomp = %d seconds' % (time()-tt00))



