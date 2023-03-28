import numpy as np
import pandas as pd
from scipy import special
from lo_tools import zfun, Lfun
from lo_user_tools import llxyfun as lxf

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

s_dict = {'THETA_S': 4, 'THETA_B': 2, 'TCLINE': 10, 'N': 30, 'VTRANSFORM': 2, 'VSTRETCHING': 4}
base_gridname = 'sillfine'
base_tag = 'v0'
gridname = 'sillfine'

# analytical model estuary with sills
# narrower, longer sill
# reduce grid domain on land side

# dch
dch = gfun.default_choices()
dch['analytical'] = True
dch['nudging_edges'] = ['north', 'south', 'west']
dch['use_z_offset'] = False
dch['z_offset'] = 0.0
dch['t_dir'] = 'BLANK'
dch['t_list'] = ['BLANK']
dch['z_land'] = 0

# fixed bathymetry parameters
L_basin = 40000 #40km long basins
W_max = 8000 # 8km wide
D_max = 200 # 200m max depth
TL_sill = 2000 # steepness of sills (transition length in m, approx 10% slope)
TL_side = 2000 # steepness of sides and end (transition length in m, approx 10% slope)

# variable bathymetry parameters
L_sill = 8000 # length of sill (length of flat section)
HR_sill = 0.75 # height ratio of sills
CR_sill = 0.5 # constriction ratio at sills

# calculate additional constants
L_estuary = (2*L_basin)+L_sill
x_sill = L_estuary/2
stretch_sill = TL_sill/4
shift_sill = L_sill/2
stretch_side = TL_side/4

# resolution constants
res_est=320 #resolution inside estuary
res_ocn=2500 #resolution outside estuary

# make grid
x_list_est_1 = np.flip(np.arange(x_sill, -res_est, -res_est)[1:])
x_list_est_2 = np.arange(x_sill, L_estuary+res_est, res_est)
x_list_est = np.concatenate((x_list_est_1, x_list_est_2), axis=None)
lon_list_est = lxf.x2lon(x_list_est,0,45)

y_list_est_1 = np.flip(np.arange(0, -12000, -res_est)[1:])
y_list_est_2 = np.arange(0, 12000, res_est)
y_list_est = np.concatenate((y_list_est_1, y_list_est_2), axis=None)
lat_list_est = lxf.y2lat(y_list_est,45)

lon_list_1 = [-lon_list_est[0], 4]
x_res_list_1 = [res_est, res_ocn]
lat_list_1 = [-lat_list_est[0], -43]
y_res_list_1 = [res_est, res_ocn]
Lon_vec_1, Lat_vec_1 = gfu.stretched_grid(lon_list_1, x_res_list_1, lat_list_1, y_res_list_1)
Lon_vec_1=-np.flip(Lon_vec_1[1:])
Lat_vec_1=-np.flip(Lat_vec_1[1:])

lon_list_2 = [lon_list_est[-1], 2]
x_res_list_2 = [res_est, res_ocn]
lat_list_2 = [lat_list_est[-1], 47]
y_res_list_2 = [res_est, res_ocn]
Lon_vec_2, Lat_vec_2 = gfu.stretched_grid(lon_list_2, x_res_list_2, lat_list_2, y_res_list_2)
Lon_vec_2=Lon_vec_2[1:]
Lat_vec_2=Lat_vec_2[1:]

Lon_vec = np.concatenate((Lon_vec_1,lon_list_est,Lon_vec_2), axis=None)
Lat_vec = np.concatenate((Lat_vec_1,lat_list_est,Lat_vec_2), axis=None)

lon, lat = np.meshgrid(Lon_vec, Lat_vec)

x, y = zfun.ll2xy(lon, lat, 0, 45)

# make bathymetry by hand
#end profile
D_end = ((D_max/2)*(special.erf((-(x-L_estuary)/stretch_side)-2)))+(D_max/2)
#constrictions
W_constrict=W_max-(((CR_sill*W_max/2)*(-special.erf(np.abs((x-x_sill)/stretch_sill)-(2+shift_sill/stretch_sill))))+(CR_sill*W_max/2))
W_bottom=W_constrict-(2*TL_side)
shift_bottom=W_bottom/2
#basin shape
z_basin=((D_end/2)*(special.erf(np.abs(y/stretch_side)-(2+shift_bottom/stretch_side))))-(D_end/2)
#sill shape
z_sill=-(((HR_sill*D_max/2)*(special.erf(np.abs((x-x_sill)/stretch_sill)-(2+shift_sill/stretch_sill))))+(D_max*(1-HR_sill/2)))
#shelf
z_shelf = x * 1e-3
z_shelf = np.where(x>=0,0,z_shelf)
#combine surfaces
z_estuary=np.maximum(z_basin,z_sill)
z_estuary=np.where(z_estuary>=-0.25,0,z_estuary)
z=np.minimum(z_estuary,z_shelf)

# create a river file
Ldir = Lfun.Lstart(gridname=base_gridname, tag=base_tag)
ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag']
Lfun.make_dir(ri_dir)
gri_fn = ri_dir / 'river_info.csv'
with open(gri_fn, 'w') as rf:
    rf.write('rname,usgs,ec,nws,ratio,depth,flow_units,temp_units\n')
    rf.write('creek_sill2,,,,1.0,5.0,m3/s,degC\n')
# and make a track for the river
track_dir = ri_dir / 'tracks'
Lfun.make_dir(track_dir)
track_fn = track_dir / 'creek_sill2.p'
track_df = pd.DataFrame()
NTR = 100
track_df['lon'] = np.linspace(0,2,NTR) # OK to go past edge of domain
track_df['lat'] = 45*np.ones(NTR)
track_df.to_pickle(track_fn)
# NOTE: tracks go from ocean to land

# check for odd size of grid and trim if needed
NR, NC = lon.shape
if np.mod(NR,2) != 0:
    print('- trimming row from grid')
    lon = lon[:-1,:]
    lat = lat[:-1,:]
    z = z[:-1,:]
    z_basin = z_basin[:-1,:]
    z_sill = z_sill[:-1,:]
    z_shelf = z_shelf[:-1,:]
    z_estuary = z_estuary[:-1,:]
    x = x[:-1,:]
    y = y[:-1,:]
if np.mod(NC,2) != 0:
    print('- trimming column from grid')
    lon = lon[:,:-1]
    lat = lat[:,:-1]
    z = z[:,:-1]
    z_basin = z_basin[:,:-1]
    z_sill = z_sill[:,:-1]
    z_shelf = z_shelf[:,:-1]
    z_estuary = z_estuary[:,:-1]
    x = x[:,:-1]
    y = y[:,:-1]