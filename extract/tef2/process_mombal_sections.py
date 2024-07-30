"""
Process tef2 extractions, giving transport vs. salinity for:
volume, salt, and other variables.

PERFORMANCE: 21 seconds for test.

To test on mac:
run process_sections.py -gtx cas7_trapsV00_meV00 -ctag c0 -0 2017.07.04 -1 2017.07.06

"""

import sys
import xarray as xr
import numpy as np
import pickle
from time import time
import pandas as pd
from scipy.stats import binned_statistic
import seawater as sw

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
in_dir2 = out_dir0 / ('extractions_stress_' + Ldir['ds0'] + '_' + Ldir['ds1'])
in_dir3 = out_dir0 / ('extractions_rpm_' + Ldir['ds0'] + '_' + Ldir['ds1'])
in_dir4 = out_dir0 / ('extractions_coriolis_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('processed_mombal_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

sect_list = [item.name for item in in_dir.glob('*.nc')]
if Ldir['testing']:
    sect_list = ['jdf3.nc']
sn_list = [item.replace('.nc','') for item in sect_list]
    
# make vn_list by inspecting the first section
ds = xr.open_dataset(in_dir / sect_list[0])
vn_list = [item for item in ds.data_vars \
    if (len(ds[item].dims) == 3) and (item not in ['vel','DZ'])]
ds.close()

# create the dict S
S_info_dict = Lfun.csv_to_dict(Ldir['grid'] / 'S_COORDINATE_INFO.csv')
S = zrfun.get_S(S_info_dict)

# get the grid file
gds = xr.open_dataset(Ldir['grid'] / 'grid.nc')
lou = gds.lon_u[0,:].values
lau = gds.lat_u[:,0].values
lov = gds.lon_v[0,:].values
lav = gds.lat_v[:,0].values
lor = gds.lon_rho.values
lar = gds.lat_rho.values
# plon, plat = pfun.get_plon_plat(lor,lar)
# hh = gds.h.values
# maskr = gds.mask_rho.values
# zm = -np.ma.masked_where(maskr==0, hh)

# get sect_df with the section point locations
sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

print('\nProcessing TEF extraction:')
print(str(in_dir))

tt00 = time()

for sn in sn_list:
    tt0 = time()
    ext_fn=sn+'.nc'
    print(ext_fn)

    # name output file
    out_fn = ext_fn

    # load fields
    ds = xr.open_dataset(in_dir / ext_fn)
    V = dict()
    for vn in vn_list:
        V[vn] = ds[vn].to_numpy()
    #V['salt2'] = V['salt']*V['salt']
    q = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['vel'].to_numpy()
    V['q'] = q
    ot = ds['time'].to_numpy()
    zeta = ds['zeta'].to_numpy()
    u_hourly = ds.vel.values #to get stress
    h = ds.h.values
    zr, zw = zrfun.get_z(h, 0*h, S) # packed (z,p)
    dz = np.diff(zw,axis=0) #does not change with time
    dzr = np.diff(zr,axis=0)
    ds.close()

    #calculate storage term
    dudt=(u_hourly[2:,:,:]-u_hourly[:-2,:,:])/3600
    dudt=np.concatenate((np.nan * np.ones((1,u_hourly.shape[1],u_hourly.shape[2])),dudt,np.nan * np.ones((1,u_hourly.shape[1],u_hourly.shape[2]))),axis=0)
    V['dudt']=dudt

    #load rpm fields
    ds3 = xr.open_dataset(in_dir3 / ext_fn)
    zetarp = ds3.zetarp.values
    zetarm = ds3.zetarm.values
    saltrp = ds3.saltrp.values
    saltrm = ds3.saltrm.values
    DZ = ds3['DZ'].values #note DZ is time varying and dz is not (calculated using zeta=0)


    # get distance array dx to use with the rpm fields (should be basically the same for each cell)
    sdf = sect_df.loc[sect_df.sn==sn,:]
    lorp = lor[sdf.jrp,sdf.irp]
    lorm = lor[sdf.jrm,sdf.irm]
    larp = lar[sdf.jrp,sdf.irp]
    larm = lar[sdf.jrm,sdf.irm]
    dxrpm = []
    for i in range(len(lorp)):
        dx, ang = sw.dist([larm[i],larp[i]],[lorm[i],lorp[i]],units='km')
        dxrpm.append(dx[0]*1000) #put into m

    # calculate non tidally averaged ds/dx and dzeta/dx
    # dsdx=zfun.lowpass((saltrp-saltrm)/dxrpm, f='godin')[pad:-pad+1:24, :]
    # dzetadx=zfun.lowpass((zetarp-zetarm)/dxrpm, f='godin')[pad:-pad+1:24, :]
    dsdx=(saltrp-saltrm)/dxrpm
    dzetadx=(zetarp-zetarm)/dxrpm

    #calculate pressure gradient
    g=9.81
    beta=7.7e-4
    pgzeta = g*dzetadx[:,np.newaxis,:]
    pgs = beta*g*np.flip(np.cumsum(np.flip(dz*dsdx,axis=1),axis=1),axis=1) #need to flip since we want to sum starting from the top which is the last element
    pgscorr = beta*g*( np.flip(np.cumsum(np.flip(dz*dsdx,axis=1),axis=1),axis=1) - 0.5*dz*dsdx)
    pg=pgzeta+pgs 
    pgcorr=pgzeta+pgscorr

    V['pg']=pgcorr

    # load stress fields
    ds2 = xr.open_dataset(in_dir2 / ext_fn)
    AKv_hourly = ds2.AKv.values
    bustr_hourly = ds2.bustr.values
    RHO0=1023.7 #mean density from BLANK.in file for model
    DZR = ds2['DZR'].values
    
    # V = dict()
    # for vn in vn_list:
    #     V[vn] = ds[vn].to_numpy()
    # #V['salt2'] = V['salt']*V['salt']
    # q = ds['dd'].to_numpy() * ds['DZ'].to_numpy() * ds['vel'].to_numpy()
    # V['q'] = q
    # ot = ds['time'].to_numpy()
    # zeta = ds['zeta'].to_numpy()
    ds2.close()

    # calculate stress section
    dudz_hourly = np.diff(u_hourly,axis=1)/DZR
    ustr_hourly = AKv_hourly*dudz_hourly
    ustr_hourly_full = np.concatenate((bustr_hourly[:,np.newaxis,:]/RHO0,ustr_hourly,np.zeros((bustr_hourly.shape[0],1,bustr_hourly.shape[1]))),axis=1)
    #ustr = zfun.lowpass(ustr_hourly_full, f='godin')[pad:-pad+1:24, :] #tidally average and subsample daily #DO NOT TIDALLY AVERAGE THIS WILL BE IN THE BULK STEP
    #dustrdz = np.diff(ustr,axis=1)/dz
    dustrdz = np.diff(ustr_hourly_full,axis=1)/DZ #since we are not tidally averaging, use ustr_hourly_full and the full time-varying DZ
    dustrdz[:,0,:]=np.nan #for now, set bottom element to nan since bustr not behaving as expected
    V['stressdiv']=dustrdz
    
    # load v section
    ds4 = xr.open_dataset(in_dir4 / ext_fn)
    v_hourly = ds4.v.values
    f = ds4.f.values
    coriolis = f[:,np.newaxis,:] * v_hourly
    V['coriolis']=coriolis


    # make arrays of property transport
    QV = dict()
    for vn in V.keys():
        if vn == 'q':
            QV[vn] = q
        else:
            QV[vn] = q*V[vn]
    NT, NZ, NX = q.shape
    
    # define salinity bins
    if Ldir['testing']:
        NS = 36 # number of salinity bins
    else:
        NS = 1000 # number of salinity bins
    S_low = 0
    S_hi = 36
    sedges = np.linspace(S_low, S_hi, NS+1)
    sbins = sedges[:-1] + np.diff(sedges)/2

    # TEF variables
    Omat = np.zeros((NT, NS))
    TEF = dict()
    for vn in QV.keys():
        TEF[vn] = Omat.copy()
    # other variables
    omat = np.zeros(NT)
    qnet = omat.copy()
    fnet = omat.copy()
    ssh = omat.copy()
    g = 9.8
    rho = 1025

    # Process into salinity bins.
    for tt in range(NT):
        sf = V['salt'][tt,:,:].squeeze().flatten()
        
        for vn in QV.keys():
            XF = QV[vn][tt,:,:].squeeze().flatten()
            if vn == 'q':
                # also keep track of volume transport
                #qnet[tt] = XF.sum() #jx
                qnet[tt] = np.nansum(XF) #jx
                # and tidal energy flux
                #zi = zeta[tt,:].squeeze()
                zi = zeta[tt,:]
                zi[np.isnan(QV[vn][tt,0,:])] = np.nan # jx, if q = NaN, zi = NaN
                #ssh[tt] = zi.mean() #jx
                ssh[tt] = np.nanmean(zi) #jx
                fnet[tt] = g * rho * ssh[tt] * qnet[tt]
                
            # scipy.stats.binned_statistic(x, values, statistic='mean', bins=10, range=None)
            XF[np.isnan(XF)] = 0 #jx
            TEF[vn][tt,:] = binned_statistic(sf, XF, statistic='sum', bins=NS, range=(S_low,S_hi)).statistic
            
            if Ldir['testing'] and (tt==10) and (vn=='salt'):
                print(TEF[vn][tt,:])
    
    TEF['qnet'] = qnet
    TEF['fnet'] = fnet
    TEF['ssh'] = ssh
    
    # Pack results in a Dataset and then save to NetCDF
    ds = xr.Dataset(coords={'time': ot,'sbins': sbins})
    for vn in ['qnet','fnet','ssh']:
        ds[vn] = (('time'), TEF[vn])
    # for vn in vn_list + ['q','salt2']: #CHANGE THIS
    for vn in vn_list + ['q','dudt','stressdiv','pg','coriolis']: #CHANGE THIS
        ds[vn] = (('time','sbins'), TEF[vn])
    # save it to NetCDF
    ds.to_netcdf(out_dir / out_fn)
    if Ldir['testing']:
        pass
    else:
        ds.close()
    
    print('  elapsed time for section = %d seconds' % (time()-tt0))
    sys.stdout.flush()
    
    

print('\nTotal elapsed time = %d seconds' % (time()-tt00))



