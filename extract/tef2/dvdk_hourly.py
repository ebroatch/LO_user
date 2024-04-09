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

from lo_tools import Lfun, zrfun, zfun
from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

# import tef_fun

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('extractions_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('dvdk_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

sect_list = [item.name for item in in_dir.glob('*.nc')]
if Ldir['testing']:
    sect_list = ['jdf3.nc']
    
# make vn_list by inspecting the first section
ds = xr.open_dataset(in_dir / sect_list[0])
vn_list = [item for item in ds.data_vars \
    if (len(ds[item].dims) == 3) and (item not in ['vel','DZ'])]
ds.close()

print('\nProcessing DvdK extraction:')
print(str(in_dir))

tt00 = time()
pad=36

for ext_fn in sect_list:
    tt0 = time()
    print(ext_fn)

    # name output file
    # out_fn = ext_fn.replace('.nc','.p')
    out_fn=ext_fn #output is now also a .nc file

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
    ds.close()
    NT, NZ, NX = q.shape

    #averaging will be different here for DvdK vs standard decomposition
    A = np.sum(dA, axis=(1,2)) #shape NT

    sbar=np.sum(sdA, axis=(1,2))/A #shape NT
    s0 = zfun.lowpass(sbar, f='godin')[pad:-pad+1] #shape NT-72
    s1=sbar[pad:-pad+1]-s0 #shape NT-72
    sprime=s[pad:-pad+1, :, :]-s1[:, np.newaxis, np.newaxis]-s0[:, np.newaxis, np.newaxis] #shape NT-72,NZ,NX newaxis for broadcasting

    ubar=np.sum(u*dA, axis=(1,2))/A
    u0 = zfun.lowpass(ubar, f='godin')[pad:-pad+1] #shape NT-72
    u1=ubar[pad:-pad+1]-u0 #shape NT-72
    uprime=u[pad:-pad+1, :, :]-u1[:, np.newaxis, np.newaxis]-u0[:, np.newaxis, np.newaxis] #shape NT-72,NZ,NX newaxis for broadcasting

    Q=np.sum(q, axis=(1,2)) #shape NT
    Q0=zfun.lowpass(Q, f='godin')[pad:-pad+1] #shape NT-72
    Q1=Q[pad:-pad+1]-Q0 #shape NT-72

    F0=Q0*s0 #shape NT-72
    F1=zfun.lowpass(Q1*s1,f='godin')[pad:-pad+1] #shape NT-144
    F2=zfun.lowpass(np.sum(uprime*sprime*dA[pad:-pad+1, :, :], axis=(1,2)), f='godin')[pad:-pad+1] #shape NT-144, NZ,NX

    F = zfun.lowpass(np.sum(u*s*dA, axis=(1,2)), f='godin')[pad:-pad+1] #shape NT-72
    ssh = np.mean(zeta, axis=1)[pad:-pad+1] #shape NT-72
    ssh_lp = zfun.lowpass(np.mean(zeta, axis=1), f='godin')[pad:-pad+1] #shape NT-72
    
    #Lateral/vertical decomposition
    uprimeL = (np.sum(uprime*dz[pad:-pad+1, :, :],axis=1))/(H[pad:-pad+1, :]) #shape NT-72,NX
    uprimeV = uprime - np.expand_dims(uprimeL,axis=1) #shape NT-72,NZ,NX expand_dims for broadcasting
    sprimeL = (np.sum(sprime*dz[pad:-pad+1, :, :],axis=1))/(H[pad:-pad+1, :]) #shape NT-72,NX
    sprimeV = sprime - np.expand_dims(sprimeL,axis=1) #shape NT-72,NZ,NX expand_dims for broadcasting

    F2L = zfun.lowpass(np.sum(uprimeL*sprimeL*H[pad:-pad+1, :]*dd, axis=1), f='godin')[pad:-pad+1] #shape NT-144
    F2V = zfun.lowpass(np.sum(uprimeV*sprimeV*dA[pad:-pad+1, :, :], axis=(1,2)), f='godin')[pad:-pad+1] #shape NT-144

    #add saving to netcdf
    #dict of variables with compatible padding removal
    dvdk = dict()
    dvdk['A']=(A[pad:-pad+1])[pad:-pad+1]
    dvdk['sbar']=(sbar[pad:-pad+1])[pad:-pad+1]
    dvdk['s0']=s0[pad:-pad+1]
    dvdk['s1']=s1[pad:-pad+1]
    dvdk['ubar']=(ubar[pad:-pad+1])[pad:-pad+1]
    dvdk['u0']=u0[pad:-pad+1]
    dvdk['u1']=u1[pad:-pad+1]
    dvdk['Q']=(Q[pad:-pad+1])[pad:-pad+1]
    dvdk['Q0']=Q0[pad:-pad+1]
    dvdk['Q1']=Q1[pad:-pad+1]
    dvdk['F0']=F0[pad:-pad+1]
    dvdk['F1']=F1
    dvdk['F2']=F2
    dvdk['F2L']=F2L
    dvdk['F2V']=F2V
    dvdk['F']=F[pad:-pad+1]
    dvdk['ssh']=ssh[pad:-pad+1]
    dvdk['ssh_lp']=ssh_lp[pad:-pad+1]
    dvdk['ot']=(ot[pad:-pad+1])[pad:-pad+1]

    #Saving as .nc dataset based on bulk_calc
    sdlist=['A','sbar','s0','s1','ubar','u0','u1','Q','Q0','Q1','F0','F1','F2','F2L','F2V','F','ssh','ssh_lp']
    # Pack results in a Dataset and then save to NetCDF
    ds = xr.Dataset(coords={'time': dvdk['ot']})
    for vn in sdlist:
        ds[vn] = (('time'), dvdk[vn])
    # save it to NetCDF
    ds.to_netcdf(out_dir / out_fn)
    if Ldir['testing']:
        pass
    else:
        ds.close()
    print('  elapsed time for section = %d seconds' % (time()-tt0))
    sys.stdout.flush()

print('\nTotal elapsed time = %d seconds' % (time()-tt00))



