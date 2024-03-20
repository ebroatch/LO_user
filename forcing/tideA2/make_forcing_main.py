"""
This is the main program for making the TIDE forcing file. TRIPLE THE AMPLITUDE FOR IDEALIZED MODEL

Test on mac in ipython:

run make_forcing_main.py -g ae0 -r backfill -d 2020.01.01 -f tideA1 -test True

"""

from pathlib import Path
from datetime import datetime, timedelta

from lo_tools import forcing_argfun2 as ffun

Ldir = ffun.intro() # this handles all the argument passing
result_dict = dict()
result_dict['start_dt'] = datetime.now()

# ****************** CASE-SPECIFIC CODE *****************

date_string = Ldir['date_string']
out_dir = Ldir['LOo'] / 'forcing' / Ldir['gridname'] / ('f' + date_string) / Ldir['frc']

import xarray as xr
import numpy as np
from lo_tools import zrfun

if Ldir['testing']:
    from importlib import reload
    reload(zrfun)

out_fn = out_dir / 'tides.nc'
out_fn.unlink(missing_ok=True)
grid_fn = Ldir['grid'] / 'grid.nc'

G = zrfun.get_basic_info(grid_fn, only_G=True)
NR, NC = G['lon_rho'].shape

# Make a dict associating a constituent with a tuple of:
# (period [hours], amplitude[m])
#cons_dict =  {'m2':(12.42, 0.75), 's2':(12, 0.25)} #original amplitude
if Ldir['gridname']=='sill5km':
    cons_dict =  {'m2':(12.42, 2.25), 's2':(12, 0.75)} #TRIPLE AMPLITUDE
elif Ldir['gridname']=='sill20km':
    cons_dict =  {'m2':(12.42, 2.06422), 's2':(12, 0.68807)} #reduce by factor of 1.09
elif Ldir['gridname']=='sill80km':
    cons_dict =  {'m2':(12.42, 1.36363), 's2':(12, 0.45454)} #reduce by factor of 1.65
NP = len(cons_dict)
mr3 = np.ones((NP, NR, NC)) * G['mask_rho'].reshape((1, NR, NC))
omat = np.zeros((NP, NR, NC))
omat[mr3==0] = np.nan

# Create the period and Eamp arrays
ds = xr.Dataset()
p_list = []
Eamp = omat.copy()
ii = 0
for c in cons_dict.keys():
    p_list.append(cons_dict[c][0])
    a = cons_dict[c][1]
    Eamp[ii,:,:] = a * np.ones((1, NR, NC))
    ii += 1
p_vec = np.array(p_list)
Eamp[mr3==0] = np.nan

# Write the period coordinate
vn = 'tide_period'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
ds[vn] = (('tide_period',), p_vec)
ds[vn].attrs['units'] = vinfo['units']
ds[vn].attrs['long_name'] = vinfo['long_name']

# Write Eamp
vn = 'tide_Eamp'
vinfo = zrfun.get_varinfo(vn, vartype='climatology')
dims = ('tide_period',) + vinfo['space_dims_tup']
ds[vn] = (dims, Eamp)
ds[vn].attrs['units'] = vinfo['units']
ds[vn].attrs['long_name'] = vinfo['long_name']

# Write all other fields
for vn in ['tide_Ephase', 'tide_Cangle', 'tide_Cphase', 'tide_Cmax', 'tide_Cmin']:
    vinfo = zrfun.get_varinfo(vn, vartype='climatology')
    dims = ('tide_period',) + vinfo['space_dims_tup']
    ds[vn] = (dims, omat.copy())
    ds[vn].attrs['units'] = vinfo['units']
    ds[vn].attrs['long_name'] = vinfo['long_name']

# Compress and save to NetCDF
Enc_dict = {vn:zrfun.enc_dict for vn in ds.data_vars}
ds.to_netcdf(out_fn, encoding=Enc_dict)
ds.close()

# -------------------------------------------------------

# test for success
if out_fn.is_file():
    result_dict['result'] = 'success' # success or fail
else:
    result_dict['result'] = 'fail'

# *******************************************************

result_dict['end_dt'] = datetime.now()
ffun.finale(Ldir, result_dict)
