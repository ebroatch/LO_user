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

t_spring = pd.Timestamp('2020-07-01 00:00:00')
t_neap = pd.Timestamp('2020-07-08 00:00:00')

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

    fig, [ax0,ax1,ax2,ax3] = plt.subplots(2, 2,figsize=(15,15))
    ax0.pcolormesh(ds['vel'].sel(time=t_spring))
    ax1.pcolormesh(ds['salt'].sel(time=t_spring))
    ax2.pcolormesh(ds['vel'].sel(time=t_neap))
    ax3.pcolormesh(ds['salt'].sel(time=t_neap))

    axs[0].set_title('u spring')
    axs[1].set_title('salt spring')
    axs[2].set_title('u neap')
    axs[3].set_title('salt neap')
    
    fig.suptitle(Ldir['sect_name'])
    plt.savefig(out_dir / (Ldir['sect_name'] + '.png'))
    ds.close()
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



