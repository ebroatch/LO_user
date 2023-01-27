"""
Generic code to plot any mooring extraction
"""
from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun

import xarray as xr
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import cycler
plt.rcParams['axes.grid'] = True
# colors = plt.cm.YlGnBu(np.linspace(0.3, 1, 5))
# plt.rcParams['axes.prop_cycle'] = cycler.cycler('color', colors)
plt.rcParams['axes.prop_cycle'] = cycler.cycler(color = ['tab:red', 'tab:orange', 'tab:green', 'tab:blue', 'tab:purple'])

Ldir = Lfun.Lstart()

# choose the file
num_fn = input("Enter number of moorings: ")
num_fn = int(num_fn)
in_dir0 = Ldir['LOo'] / 'extract'
gtagex = Lfun.choose_item(in_dir0, tag='', exclude_tag='', itext='** Choose gtagex from list **')
in_dir = in_dir0 / gtagex / 'moor'
# you can choose either a file or a directory
moor_name = Lfun.choose_item(in_dir, itext='** Choose mooring extraction or folder from list **')
moor_item = in_dir / moor_name
if moor_item.is_file() and moor_name[-3:]=='.nc':
    moor_fn = moor_item
    num_fn = 1
    moor_label = input("Enter label for mooring")

fig, axs = plt.subplots(4,1,sharex=True)
lp = True
if lp == True:
    axs[0].set_title('Tidally averaged mooring extractions')
else:
    axs[0].set_title('Mooring extractions')


for i in range(num_fn):
    if moor_item.is_dir():
        moor_name = Lfun.choose_item(moor_item, tag='.nc', itext='** Choose mooring extraction #'+str(i+1)+' from list **')
        moor_fn = moor_item / moor_name
        moor_label = input("Enter label for mooring: ")

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

    # populate lists of variables to plot
    vn2_list = ['zeta','ubar','vbar']
    vn3_list = ['salt']

    # drop missing variables
    vn2_list = [item for item in vn2_list if item in ds.data_vars]
    vn3_list = [item for item in vn3_list if item in ds.data_vars]

    # plot time series using a pandas DataFrame
    df = pd.DataFrame(index=ot)
    for vn in vn2_list:
        if lp == True:
            df[vn] = zfun.lowpass(ds[vn].values, f='godin')
        else:
            df[vn] = ds[vn].values
        
    for vn in vn3_list:
        if lp == True:
            df[vn] = zfun.lowpass(ds[vn][:, -1].values, f='godin')
        else:
            df[vn] = ds[vn][:, -1]

    axs[0].plot(df['zeta'], label=moor_label)
    axs[0].set_ylabel(r'$\zeta$ [m]')
    axs[0].text(.02, .9, 'Sea surface height', c='k', weight='bold', transform=axs[0].transAxes)

    axs[1].plot(df['salt'], label=moor_label)
    axs[1].set_ylim(bottom=0,top=30)
    axs[1].set_ylabel(r'S [psu]')
    axs[1].text(.02, .9, 'Surface salinity', c='k', weight='bold', transform=axs[1].transAxes)

    axs[2].plot(df['ubar'], label=moor_label) 
    axs[2].set_ylabel(r'u [m/s]')
    axs[2].text(.02, .9, 'ubar', c='k', weight='bold', transform=axs[2].transAxes)

    axs[3].plot(df['vbar'], label=moor_label)
    axs[3].set_ylabel(r'$\bar{v}$ [m/s]')
    axs[3].text(.02, .9, 'vbar', c='k', weight='bold', transform=axs[3].transAxes)

    if lp == True:
        axs[0].set_ylim(bottom=0, top=0.5)
        axs[2].set_ylim(bottom=-0.1, top=0.4)
        axs[3].set_ylim(bottom=-0.08, top=0.08)
    else:
        axs[0].set_ylim(bottom=-2, top=3)
        axs[2].set_ylim(bottom=-0.8, top=1.4)
        axs[3].set_ylim(bottom=-0.2, top=0.2)


    # ax = fig.add_subplot(412)
    # cs = ax.pcolormesh(TT, ZZ, TH, cmap=cm.thermal, vmin=5, vmax=15)
    # fig.colorbar(cs)
    # ax.set_ylim(top=5)
    # ax.set_ylabel('Z [m]')
    # ax.text(.05, .1, 'Potential Temperature [deg C]', c='w', weight='bold', transform=ax.transAxes)
    # #ax.set_xlabel('Time [days from start of record]')

    # ax = fig.add_subplot(413)
    # cs = ax.pcolormesh(TT, ZZ, U, cmap=cm.balance, vmin=-4, vmax=4)
    # fig.colorbar(cs)
    # ax.set_ylim(top=5)
    # ax.set_ylabel('Z [m]')
    # ax.text(.05, .1, 'u [m/s]', c='w', weight='bold', transform=ax.transAxes)
    # #ax.set_xlabel('Time [days from start of record]')

    # ax = fig.add_subplot(414)
    # cs = ax.pcolormesh(TT, ZZ, V, cmap=cm.balance, vmin=-1, vmax=1)
    # fig.colorbar(cs)
    # ax.set_ylim(top=5)
    # ax.set_ylabel('Z [m]')
    # ax.text(.05, .1, 'v [m/s]', c='w', weight='bold', transform=ax.transAxes)
    # ax.set_xlabel('Time [days from start of record]')
axs[0].legend(loc='best')
# axs[1].legend(loc='best')
# axs[2].legend(loc='best')
# axs[3].legend(loc='best')
axs[0].set_xlim(left=pd.Timestamp('2020-01-01'), right=pd.Timestamp('2020-04-01'))
plt.show()