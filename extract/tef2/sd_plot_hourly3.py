"""
Plot bulk fluxes as a time series.

To test on mac:
run bulk_plot -gtx cas6_v00_uu0m -ctag c0 -0 2022.01.01 -1 2022.12.31 -test True


"""
import sys
import matplotlib.pyplot as plt
import numpy as np
import pickle
from time import time
import pandas as pd
import xarray as xr

from lo_tools import Lfun, zfun
from lo_tools import plotting_functions as pfun
import flux_fun

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('standard_decomp_' + Ldir['ds0'] + '_' + Ldir['ds1'])
in_dir2 = out_dir0 / ('processed_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('sd_plots_hourly3_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

# sect_list = [item.name for item in in_dir.glob('*.p')]
# if Ldir['testing']:
#     sect_list = ['ss2.p']
# sect_list = ['a1.p','a3.p','b1.p','b3.p','b5.p','c3.p']
# sect_label = ['a1','a3','b1','b3','b5','c3']

# sect_list = ['b1.p','b2.p','b3.p','b4.p','b5.p']
# sect_nclist = ['b1.nc','b2.nc','b3.nc','b4.nc','b5.nc']
# sect_label = ['b1','b2','b3','b4','b5']
# plot_color = ['tab:red','tab:orange','tab:green','tab:cyan','tab:blue','tab:purple']

sect_list = ['a1.p','a2.p','a3.p','a4.p','a5.p','b1.p','b2.p','b3.p','b4.p','b5.p','c1.p','c2.p','c3.p','c4.p','c5.p']
sect_label = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
plot_color = ['pink','tab:pink','hotpink','deeppink','magenta','tab:red','tab:orange','tab:green','tab:cyan','tab:blue','mediumslateblue','tab:purple','blueviolet','rebeccapurple','indigo']

#sect_list = ['a1.p','a3.p','b1.p','b2.p','b3.p','b4.p','b5.p','c3.p']
#sect_label = ['a1','a3','b1','b2','b3','b4','b5','c3']
#plot_color = ['tab:red','tab:orange','tab:olive','tab:green','tab:cyan','tab:blue','tab:purple','tab:pink']
#plot_color = ['k','tab:gray','tab:red','tab:orange','tab:green','tab:cyan','tab:blue','tab:brown']


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

# PLOTTING
#plot_color = ['lightblue','tab:cyan','dodgerblue','tab:blue','blue','gold','goldenrod','xkcd:yellow orange','tab:orange','peru','pink','tab:pink','mediumvioletred','tab:red','maroon']
#plot_label = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
# p_color = ['tab:blue','tab:orange','tab:red']
# m_color = ['tab:cyan','xkcd:yellow orange','tab:pink']
# label_in = ['a3 in','b3 in','c3 in']
# label_out = ['a3 out','b3 out','c3 out']
fs = 12
plt.close('all')
pfun.start_plot(fs=fs, figsize=(21,10))

fig, axs = plt.subplots(2, 1, sharex=True,figsize=(15,7.7),gridspec_kw={'height_ratios': [4,1]})
# fig = plt.figure()   
# ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
# ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
# ax3 = plt.subplot2grid((1,3), (0,2)) # map

for i in range(len(sect_list)):
    sect_name = sect_list[i]
    SD = pickle.load(open(in_dir / sect_name, 'rb'))

    FR = SD['FR']
    FE = SD['FE']
    FT = SD['FT']
    FTL = SD['FTL']
    FTV = SD['FTV']
    F = SD['F']
    ot = SD['ot']
                    
    # labels and colors
    # ylab_dict = {'Q': r'Transport $[10^{3}\ m^{3}s^{-1}]$',
    #             'salt': r'Salinity $[g\ kg^{-1}]$'}
    # ylab_dict = {'Q': r'$Q_{in}\ [10^{3}\ m^{3}s^{-1}]$',
    #             'salt': r'$s_{in}\ [g\ kg^{-1}]$',
    #             'deltas': r'$\Delta s\ [g\ kg^{-1}]$'}
    # p_color = 'r'
    # m_color = 'b'
    lw = 2
        
    # fig = plt.figure()
    
    # ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
    # ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
    # ax3 = plt.subplot2grid((1,3), (0,2)) # map
    
    #ot = bulk['ot'] # (same as tef_df.index)
    
    # axs[0].plot(ot,FR, color='tab:blue', linewidth=lw, label=sect_label[i])
    # axs[2].plot(ot,FE, color='tab:green', linewidth=lw, label=sect_label[i])
    # axs[3].plot(ot,FT, color='tab:red', linewidth=lw, label=sect_label[i])
    axs[1].plot(ot, FR, color=plot_color[i], linewidth=lw, label=sect_label[i])
    axs[0].plot(ot, FE+FT, color=plot_color[i], linewidth=lw, label=sect_label[i]+' FE+FT')
    axs[0].plot(ot, FE, linestyle = '--', color=plot_color[i], linewidth=lw, label=sect_label[i]+' FE')
    axs[0].plot(ot, FT, linestyle = ':', color=plot_color[i], linewidth=lw, label=sect_label[i]+' FT')
    # axs[i].plot(ot,FTL, linestyle = '--', color='tab:pink', linewidth=lw, label=r'$F_{TL}$')
    # axs[i].plot(ot,FTV, linestyle = ':', color='tab:orange', linewidth=lw, label=r'$F_{TV}$')

axs[0].grid(True)
axs[1].grid(True) 
#axs[2].grid(True)   
axs[1].set_ylabel(r'$FR [m^{3}s^{-1} g\ kg^{-1}]$')
axs[0].set_ylabel(r'$Salt\ flux [m^{3}s^{-1} g\ kg^{-1}]$')
#axs[2].set_ylabel(r'$F_{T}[m^{3}s^{-1} g\ kg^{-1}]$')
axs[1].set_ylim(-3.2e4,-2.8e4)
axs[0].set_ylim(-2e4,5e4)
#axs[2].set_ylim(-2e4,4e4)
axs[0].legend(loc='lower right')
axs[2].set_xlim(pd.Timestamp('2020-06-01'), pd.Timestamp('2020-07-31'))
#plt.suptitle('Standard decomposition')
plt.savefig(out_dir / ('sd_plot_hourly3.png'))
plt.close()

#scatter plot
fig, ax = plt.subplots(1, 1, figsize=(10,10))
for i in range(len(sect_list)):
    sect_name = sect_list[i]
    #sect_ncname = sect_nclist[i]
    SD = pickle.load(open(in_dir / sect_name, 'rb'))

    FR = SD['FR']
    FE = SD['FE']
    FT = SD['FT']
    FTL = SD['FTL']
    FTV = SD['FTV']
    F = SD['F']
    ot = SD['ot']

    ds2= pickle.load(open(in_dir2 / sect_name, 'rb'))
    #ds2 = xr.open_dataset(in_dir2 / sect_name) #change in future
    #qprism=zfun.lowpass(np.abs(ds2['qnet'].values-zfun.lowpass(ds2['qnet'].values, f='godin')), f='godin')/2
    qprism=zfun.lowpass(np.abs(ds2['qnet']-zfun.lowpass(ds2['qnet'], f='godin')), f='godin')/2
    pad=36
    qprism=(qprism[pad:-pad+1])[pad:-pad+1]

    ax.plot(qprism, FT/(FE+FT), color=plot_color[i], linewidth=lw, label=sect_label[i])

ax.grid(True)
ax.set_ylabel('FT/(FE+FT)')
ax.set_xlabel('Qprism')
ax.legend(loc='lower right')

plt.savefig(out_dir / ('sd_plot_scatter.png'))
plt.close()

