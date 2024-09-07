"""
Plot bulk fluxes as a time series.

To test on mac:
run bulk_plot -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31 -test True


"""
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pickle
from time import time
import pandas as pd
import xarray as xr

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import flux_fun
import tef_fun

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('standard_decomp_' + Ldir['ds0'] + '_' + Ldir['ds1'])
in_dir2 = out_dir0 / ('processed_' + Ldir['ds0'] + '_' + Ldir['ds1'])
in_dir3 = out_dir0 / ('bulk_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('sd_plots_fine_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

sect_list = [item.name for item in in_dir.glob('*.nc')]
# if Ldir['testing']:
#     sect_list = ['ss2.p']
# sect_list = ['a1.p','a3.p','b1.p','b3.p','b5.p','c3.p']
# sect_label = ['a1','a3','b1','b3','b5','c3']
# sect_list = ['b1.p','b2.p','b3.p','b4.p','b5.p']
# sect_nclist = ['b1.nc','b2.nc','b3.nc','b4.nc','b5.nc']
# sect_label = ['b1','b2','b3','b4','b5']

# plot_color = ['tab:red','tab:orange','tab:green','tab:cyan','tab:blue','tab:purple']
# plot_color = ['tab:cyan',plt.cm.Dark2(0),'tab:olive',plt.cm.Dark2(5),'tab:brown']
# plot_color_light = [plt.cm.tab20(19),plt.cm.Set2(0),plt.cm.tab20(17),plt.cm.Dark2(5),plt.cm.tab20(11)]

# sect_list = ['a1.p','a3.p','b1.p','b2.p','b3.p','b4.p','b5.p','c3.p']
# sect_label = ['a1','a3','b1','b2','b3','b4','b5','c3']
# #plot_color = ['tab:red','tab:orange','tab:olive','tab:green','tab:cyan','tab:blue','tab:purple','tab:pink']
# plot_color = ['k','tab:gray','tab:red','tab:orange','tab:green','tab:cyan','tab:blue','tab:brown']
#sect_list = ['a1.p','a2.p','a3.p','a4.p','a5.p','b1.p','b2.p','b3.p','b4.p','b5.p','c1.p','c2.p','c3.p','c4.p','c5.p']

# grid info
g = xr.open_dataset(Ldir['grid'] / 'grid.nc')
h = g.h.values
h[g.mask_rho.values==0] = np.nan
xrho = g.lon_rho.values
yrho = g.lat_rho.values
xp, yp = pfun.get_plon_plat(xrho,yrho)
xu = g.lon_u.values
yu = g.lat_u.values
xv = g.lon_v.values
yv = g.lat_v.values

# get sect_df with the section point locations
sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

#open first section to get dimensions
sect_name = 'd0'
sect_ncname = sect_name+'.nc'
SD = xr.open_dataset(in_dir / sect_ncname)
FE = SD['FE']
FT = SD['FT']
ot = SD['time']

SDfull = xr.Dataset(
    {
        'FE': (['section','time'], np.zeros((len(sect_list),len(FE)))),
        'FT': (['section','time'], np.zeros((len(sect_list),len(FT)))),
        'lon': (['section'], np.zeros(len(sect_list))),
        'xkm': (['section'], np.zeros(len(sect_list)))
    },
    coords={
        'section': np.arange(len(sect_list)) ,
        'time': ot
    },
)

for i in range(len(sect_list)):
    sect_name = 'd'+str(i)
    sect_ncname = sect_name+'.nc'

    SD = xr.open_dataset(in_dir / sect_ncname)
    #SD = pickle.load(open(in_dir / sect_name, 'rb'))
    # FR = SD['FR']
    FE = SD['FE']
    FT = SD['FT']
    # FTL = SD['FTL']
    # FTV = SD['FTV']
    # F = SD['F']
    # ot = SD['time']

    SDfull['FE'].loc[dict(section=i)]=FE
    SDfull['FT'].loc[dict(section=i)]=FT

    sdf = sect_df.loc[sect_df.sn==sect_name,:]
                    


# # PLOTTING
# #plot_color = ['lightblue','tab:cyan','dodgerblue','tab:blue','blue','gold','goldenrod','xkcd:yellow orange','tab:orange','peru','pink','tab:pink','mediumvioletred','tab:red','maroon']
# #plot_label = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
# # p_color = ['tab:blue','tab:orange','tab:red']
# # m_color = ['tab:cyan','xkcd:yellow orange','tab:pink']
# # label_in = ['a3 in','b3 in','c3 in']
# # label_out = ['a3 out','b3 out','c3 out']
# # fs = 12
# plt.close('all')
# #pfun.start_plot(fs=fs, figsize=(21,10))
# #pfun.start_plot(fs=fs, figsize=(7,10)) #narrower plot with one month xlim
# pfun.start_plot(fs=14)
# #fig, axs = plt.subplots(4, 1, sharex=True,figsize=(10,10),gridspec_kw={'height_ratios': [5,5,1,1]})
# fig, axs = plt.subplots(4, 1, sharex=True,figsize=(8,8),gridspec_kw={'height_ratios': [1,6,6,1]})
# # fig = plt.figure()   
# # ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
# # ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
# # ax3 = plt.subplot2grid((1,3), (0,2)) # map

# # labels and colors
# # ylab_dict = {'Q': r'Transport $[10^{3}\ m^{3}s^{-1}]$',
# #             'salt': r'Salinity $[g\ kg^{-1}]$'}
# # ylab_dict = {'Q': r'$Q_{in}\ [10^{3}\ m^{3}s^{-1}]$',
# #             'salt': r'$s_{in}\ [g\ kg^{-1}]$',
# #             'deltas': r'$\Delta s\ [g\ kg^{-1}]$'}
# # p_color = 'r'
# # m_color = 'b'
# lw = 2
    
# # fig = plt.figure()

# # ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
# # ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
# # ax3 = plt.subplot2grid((1,3), (0,2)) # map

# #ot = bulk['ot'] # (same as tef_df.index)

# # axs[0].plot(ot,FR, color='tab:blue', linewidth=lw, label=sect_label[i])
# # axs[2].plot(ot,FE, color='tab:green', linewidth=lw, label=sect_label[i])
# # axs[3].plot(ot,FT, color='tab:red', linewidth=lw, label=sect_label[i])
# axs[3].plot(ot, FR, color=plot_color[i], linewidth=lw, label=sect_label[i])
# axs[1].plot(ot, FE, color=plot_color[i], linewidth=lw, label=sect_label[i])
# if i==0:
#     axs[2].plot(ot, FT, color='k', linewidth=lw, label=r'$F_{T}$')
#     # axs[2].plot(ot, FTV, color='k', linewidth=lw, ls=':', label=r'$F_{TV}$')
#     axs[2].plot(ot, FTL, color='k', linewidth=lw, ls=':', label=r'$F_{TL}$')
#     axs[2].legend(loc='lower right')
# axs[2].plot(ot, FT, color=plot_color[i], linewidth=lw, label=sect_label[i])
# # axs[2].plot(ot, FTV, color=plot_color_light[i], linewidth=lw)#, ls=':')# label=sect_label[i])
# axs[2].plot(ot, FTL, color=plot_color[i], linewidth=lw, ls=':')# label=sect_label[i])
# # axs[i].plot(ot,FTL, linestyle = '--', color='tab:pink', linewidth=lw, label=r'$F_{TL}$')
# # axs[i].plot(ot,FTV, linestyle = ':', color='tab:orange', linewidth=lw, label=r'$F_{TV}$')

# axs[0].grid(True)
# axs[1].grid(True) 
# axs[2].grid(True)
# axs[3].grid(True)

# axs[0].text(0.02,0.8,'A',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[0].transAxes)
# axs[1].text(0.02,0.95,'B',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[1].transAxes)
# axs[2].text(0.02,0.95,'C',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[2].transAxes)
# axs[3].text(0.02,0.8,'D',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[3].transAxes)

# axs[1].legend(loc='upper right')

# # axs[2].set_xlim(pd.Timestamp('2020-10-01'), pd.Timestamp('2020-11-15'))
# # axs[3].set_xticks(ticks=[pd.Timestamp('2020-10-01'), pd.Timestamp('2020-10-15'), pd.Timestamp('2020-11-01'), pd.Timestamp('2020-11-15')])
# #axs[3].set_xlim(pd.Timestamp('2020-09-01'), pd.Timestamp('2020-12-31'))
# axs[3].set_xlim(pd.Timestamp('2020-10-01'), pd.Timestamp('2020-10-31'))
# #plt.xticks(rotation=90)

# axs[3].set_ylabel('$F_{R}$\n$[m^{3}s^{-1} g\ kg^{-1}]$')
# axs[1].set_ylabel('$F_{E}$\n$[m^{3}s^{-1} g\ kg^{-1}]$')
# axs[2].set_ylabel('$F_{T}$\n$[m^{3}s^{-1} g\ kg^{-1}]$')
# # axs[0].set_ylabel('$Q_{prism}$\n$[10^{3}\ m^{3}s^{-1}]$')
# axs[0].set_ylabel('$Q_{prism}$ (b3)\n$[10^{3}\ m^{3}s^{-1}]$')

# axs[3].xaxis.set_major_formatter(mdates.DateFormatter('%j')) #yearday
# axs[3].set_xlabel('Yearday')

# # if Ldir['gridname']=='sill20kmdeep':
# #     axs[0].set_title('20km sill')
# # elif Ldir['gridname']=='sill5km':
# #     axs[0].set_title('5km sill')
# # elif Ldir['gridname']=='sill80km':
# #     axs[0].set_title('80km sill')
# #plt.suptitle('Standard decomposition')
# # axs[0].set_title(Ldir['gtagex'])
# plt.tight_layout()

# plt.savefig(out_dir / ('sd_plot_hourly6.png'))
# plt.close()
# pfun.end_plot()
