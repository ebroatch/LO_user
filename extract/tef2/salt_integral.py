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
c_dir = Ldir['LOo'] / 'extract' / 'tef2' / ('sections' + '_' + gctag) 
out_dir = out_dir0 / ('salt_integral_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

sect_list = [item.name for item in in_dir.glob('*.nc')]
if Ldir['testing']:
    sect_list = ['a1.nc','b3.nc']

fn_list = Lfun.get_fn_list('hourly', Ldir, Ldir['ds0'], Ldir['ds1'])
N = len(fn_list)
    
# make vn_list by inspecting the first section
# ds = xr.open_dataset(in_dir / sect_list[0])
# vn_list = [item for item in ds.data_vars \
#     if (len(ds[item].dims) == 3) and (item not in ['vel','DZ'])]
# ds.close()

print('\nCalculating salt integrals:')
print(str(in_dir))

tt00 = time()
pad=36

for ext_fn in sect_list:
    tt0 = time()
    print(ext_fn)

    # name output file
    out_fn = ext_fn.replace('.nc','.p')
    sect_key = ext_fn.replace('.nc','')
    one_section_df = pd.read_pickle(c_dir / out_fn)
    sect_lon = one_section_df.loc[0,'x']

    ot_list=[]
    V_list=[]
    s_int_list=[]
    s_bar_list=[]

    # loop over history files
    for ii in range(N):
        fn = fn_list[ii]

        G, S, T = zrfun.get_basic_info(fn)
        h = G['h']
        DA = G['DX'] * G['DY']
        DA3 = DA.reshape((1,G['M'],G['L']))

        ds = xr.open_dataset(fn)
        ot = ds['ocean_time'][0].to_numpy()
        h = ds['h'][:].to_numpy()
        zeta = ds['zeta'][0,:,:].to_numpy()
        
        z_w = zrfun.get_z(h, zeta, S, only_w=True)
        dz = np.diff(z_w, axis=0)
        DV = dz * DA3

        s = ds['salt']
        #sdv = s*DV
        #dv = sdv/s #make DataArray
        dv = xr.ones_like(s)*DV
        sdv = s*dv
        v_int = dv.where(dv.lon_rho>sect_lon).sum()
        s_int = sdv.where(dv.lon_rho>sect_lon).sum()
        s_avg = s_int/v_int

        ot_list.append(ot)
        V_list.append(v_int)
        s_int_list.append(s_int)
        s_bar_list.append(s_avg)

        if ii%25==0:
            print(ii)

        ds.close()

    # SI = dict()
    # SI['ot']=np.asarray(ot_list)
    # SI['V']=np.asarray(V_list)
    # SI['s_int']=np.asarray(s_int_list)
    # SI['s_bar']=np.asarray(s_bar_list)
    # pickle.dump(SI, open(out_dir / out_fn, 'wb'))

    si_ds=xr.Dataset(data_vars = dict(V=('time', np.asarray(V_list)), s_int=('time', np.asarray(s_int_list)), s_bar=('time', np.asarray(s_bar_list))), coords = dict(time=np.asarray(ot_list)))
    si_ds.to_netcdf(out_dir / out_fn)
    si_ds.close()

    print('  elapsed time for section = %d seconds' % (time()-tt0))
    sys.stdout.flush()





    # # load fields
    # ds = xr.open_dataset(in_dir / ext_fn)
    # # V = dict()
    # # for vn in vn_list:
    # #     V[vn] = ds[vn].to_numpy()
    # # V['salt2'] = V['salt']*V['salt']

    # ot = ds['time'].to_numpy() #shape NT
    # zeta = ds['zeta'].to_numpy() #shape NT,NX
    # h = ds['h'].to_numpy() #shape NX
    # dd = ds['dd'].to_numpy() #shape NX
    # dz = ds['DZ'].to_numpy() #shape NT,NZ,NX
    # u = ds['vel'].to_numpy() #shape NT,NZ,NX
    # s = ds['salt'].to_numpy() #shape NT,NZ,NX
    
    # dA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() #shape NT,NZ,NX
    # H = np.sum(ds['DZ'].to_numpy(),axis=1) #shape NT,NX
    # q = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['vel'].to_numpy() #shape NT,NZ,NX
    # sdA = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['salt'].to_numpy() #shape NT,NZ,NX
    # # V['q'] = q
    # ds.close()
    # NT, NZ, NX = q.shape

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
    # SD['ssh_lp']=ssh[pad:-pad+1]
    # SD['ssh_lp']=ssh_lp[pad:-pad+1]
    # SD['ot']=(ot[pad:-pad+1])[pad:-pad+1]
    # pickle.dump(SD, open(out_dir / out_fn, 'wb'))
    # print('  elapsed time for section = %d seconds' % (time()-tt0))
    # sys.stdout.flush()

    # # # make arrays of property transport
    # # QV = dict()
    # # for vn in V.keys():
    # #     if vn == 'q':
    # #         QV[vn] = q
    # #     else:
    # #         QV[vn] = q*V[vn]
    
    
    # # # define salinity bins
    # # if Ldir['testing']:
    # #     NS = 36 # number of salinity bins
    # # else:
    # #     NS = 1000 # number of salinity bins
    # # S_low = 0
    # # S_hi = 36
    # # sedges = np.linspace(S_low, S_hi, NS+1)
    # # sbins = sedges[:-1] + np.diff(sedges)/2

    # # # TEF variables
    # # Omat = np.zeros((NT, NS))
    # # TEF = dict()
    # # for vn in QV.keys():
    # #     TEF[vn] = Omat.copy()
    # # # other variables
    # # omat = np.zeros(NT)
    # # qnet = omat.copy()
    # # fnet = omat.copy()
    # # ssh = omat.copy()
    # # g = 9.8
    # # rho = 1025

    # # # Process into salinity bins.
    # # # NOTE: this seems like a lot of nested loops (section/time/variable/digitization indices)
    # # # but it runs reasonably fast. Also we only have to run it once.
    # # for tt in range(NT):
            
    # #     if False:
    # #         sf = V['salt'][tt,:,:].squeeze().flatten()
    # #         # sort into salinity bins
    # #         inds = np.digitize(sf, sedges, right=True)
    # #         indsf = inds.copy().flatten()
            
    # #         for vn in QV.keys():
    # #             XI = QV[vn][tt,:,:].squeeze()
    # #             XF = zfun.fillit(XI).flatten()
    # #             XF = XF[~np.isnan(XF)]
    # #             if vn == 'q':
    # #                 # also keep track of volume transport
    # #                 qnet[tt] = XF.sum()
    # #                 # and tidal energy flux
    # #                 zi = zeta[tt,:].squeeze()
    # #                 ssh[tt] = zi.mean()
    # #                 fnet[tt] = g * rho * ssh[tt] * qnet[tt]
    # #             counter = 0
    # #             for ii in indsf:
    # #                 TEF[vn][tt, ii-1] += XF[counter]
    # #                 counter += 1
                
                
    # #             if Ldir['testing'] and (tt==10) and (vn=='salt'):
    # #                 print(TEF[vn][tt,:])
                
    # #     else:
    # #         # alternate version
            
    # #         sf = V['salt'][tt,:,:].squeeze().flatten()
            
    # #         for vn in QV.keys():
    # #             XF = QV[vn][tt,:,:].squeeze().flatten()
    # #             if vn == 'q':
    # #                 # also keep track of volume transport
    # #                 qnet[tt] = XF.sum()
    # #                 # and tidal energy flux
    # #                 zi = zeta[tt,:].squeeze()
    # #                 ssh[tt] = zi.mean()
    # #                 fnet[tt] = g * rho * ssh[tt] * qnet[tt]
                    
    # #             # scipy.stats.binned_statistic(x, values, statistic='mean', bins=10, range=None)
    # #             TEF[vn][tt,:] = binned_statistic(sf, XF, statistic='sum', bins=NS, range=(S_low,S_hi)).statistic
                
    # #             # results are identical
                
                
    # #             if Ldir['testing'] and (tt==10) and (vn=='salt'):
    # #                 print(TEF[vn][tt,:])
            
        
    
    # # TEF['ot'] = ot
    # # TEF['sbins'] = sbins
    # # TEF['qnet'] = qnet
    # # TEF['fnet'] = fnet
    # # TEF['ssh'] = ssh

print('\nTotal elapsed time = %d seconds' % (time()-tt00))



