"""
Generic code to plot any mooring extraction, using pcolor.
"""

from lo_tools import Lfun
from lo_tools import plotting_functions as pfun

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from cmocean import cm

Ldir = Lfun.Lstart()

# choose the file
in_dir0 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir0, tag='', exclude_tag='', itext='** Choose gtagex from list **')
in_dir = in_dir0 / gtagex / 'moor'
# you can choose either a file or a directory
moor_name = Lfun.choose_item(in_dir, itext='** Choose mooring extraction or folder from list **')
moor_item = in_dir / moor_name
if moor_item.is_file() and moor_name[-3:]=='.nc':
    moor_fn = moor_item
elif moor_item.is_dir():
    moor_name = Lfun.choose_item(moor_item, tag='.nc', itext='** Choose mooring extraction from list **')
    moor_fn = moor_item / moor_name

# load everything using xarray
ds = xr.open_dataset(moor_fn)
ot = ds.ocean_time.values
ot_dt = pd.to_datetime(ot)
t = (ot_dt - ot_dt[0]).total_seconds().to_numpy()
T = t/86400 # time in days from start
print('time step of mooring'.center(60,'-'))
print(t[1])
print('time limits'.center(60,'-'))
print('start ' + str(ot_dt[0]))
print('end   ' + str(ot_dt[-1]))
print('info'.center(60,'-'))
VN_list = []
for vn in ds.data_vars:
    print('%s %s' % (vn, ds[vn].shape))
    VN_list.append(vn)
    
s = ds['salt'].values
th = ds['temp'].values
u = ds['u'].values
v = ds['v'].values
#ox = ds['oxygen'].values
z = ds['z_w'].values
NT, NZ = z.shape
Z = z.mean(axis=0)
Z = Z - Z[-1] # adjust top to zero

# coordinate arrays for plotting
TT = T.reshape((NT,1))*np.ones((1,NZ))
ZZ = Z.reshape((1,NZ))*np.ones((NT,1))

# make variables at middle times, for pcolormesh
S = (s[1:,:] + s[:-1,:])/2
TH = (th[1:,:] + th[:-1,:])/2
U = (u[1:,:] + u[:-1,:])/2
V = (v[1:,:] + v[:-1,:])/2
# OX = (ox[1:,:] + ox[:-1,:])/2
# OX = OX*32/1000 # convert uM to mg/L

plt.close('all')
pfun.start_plot()
fig = plt.figure()

ax = fig.add_subplot(311)
cs = ax.pcolormesh(TT, ZZ, S, cmap=cm.haline, vmin=24, vmax=30)
fig.colorbar(cs)
ax.set_ylim(top=5)
ax.set_ylabel('Z [m]')
ax.set_title(moor_name)
ax.text(.05, .1, 'Salinity [g/kg]', c='k', weight='bold', transform=ax.transAxes)

# ax = fig.add_subplot(412)
# cs = ax.pcolormesh(TT, ZZ, TH, cmap=cm.thermal, vmin=5, vmax=15)
# fig.colorbar(cs)
# ax.set_ylim(top=5)
# ax.set_ylabel('Z [m]')
# ax.text(.05, .1, 'Potential Temperature [deg C]', c='w', weight='bold', transform=ax.transAxes)
# #ax.set_xlabel('Time [days from start of record]')

ax = fig.add_subplot(312)
cs = ax.pcolormesh(TT, ZZ, U, cmap=cm.balance, vmin=-0.5, vmax=0.5)
fig.colorbar(cs)
ax.set_ylim(top=5)
ax.set_ylabel('Z [m]')
ax.text(.05, .1, 'u [m/s]', c='w', weight='bold', transform=ax.transAxes)
#ax.set_xlabel('Time [days from start of record]')

ax = fig.add_subplot(313)
cs = ax.pcolormesh(TT, ZZ, V, cmap=cm.balance, vmin=-0.1, vmax=0.1)
fig.colorbar(cs)
ax.set_ylim(top=5)
ax.set_ylabel('Z [m]')
ax.text(.05, .1, 'v [m/s]', c='w', weight='bold', transform=ax.transAxes)
ax.set_xlabel('Time [days from start of record]')

# ax = fig.add_subplot(313)
# #cs = ax.pcolormesh(TT, ZZ, OX, cmap=cm.oxy, vmin=0, vmax=10)
# cs = ax.pcolormesh(TT, ZZ, OX, cmap='jet', vmin=0, vmax=10)
# fig.colorbar(cs)
# ax.set_ylim(top=5)
# ax.set_xlabel('Time [days from start of record]')
# ax.set_ylabel('Z [m]')
# ax.text(.05, .1, 'Dissolved Oxygen [mg/L]', c='k', weight='bold', transform=ax.transAxes)

plt.show()

pfun.end_plot()