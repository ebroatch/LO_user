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
out_dir = out_dir0 / ('sd_plots_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

# sect_list = [item.name for item in in_dir.glob('*.p')]
# if Ldir['testing']:
#     sect_list = ['ss2.p']
sect_list = ['a3.p','b3.p','c3.p']
#sect_list = ['a1.p','a2.p','a3.p','a4.p','a5.p','b1.p','b2.p','b3.p','b4.p','b5.p','c1.p','c2.p','c3.p','c4.p','c5.p']

# grid info
g = xr.open_dataset(Ldir['grid'] / 'grid.nc')
h = g.h.values
h[g.mask_rho.values==0] = np.nan
xr = g.lon_rho.values
yr = g.lat_rho.values
xp, yp = pfun.get_plon_plat(xr,yr)
xu = g.lon_u.values
yu = g.lat_u.values
xv = g.lon_v.values
yv = g.lat_v.values

# PLOTTING
#plot_color = ['lightblue','tab:cyan','dodgerblue','tab:blue','blue','gold','goldenrod','xkcd:yellow orange','tab:orange','peru','pink','tab:pink','mediumvioletred','tab:red','maroon']
#plot_label = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
plot_title = ['a3 outer','b3 sill','c3 inner']
# p_color = ['tab:blue','tab:orange','tab:red']
# m_color = ['tab:cyan','xkcd:yellow orange','tab:pink']
# label_in = ['a3 in','b3 in','c3 in']
# label_out = ['a3 out','b3 out','c3 out']
fs = 12
plt.close('all')
pfun.start_plot(fs=fs, figsize=(21,10))

fig, axs = plt.subplots(3, 1, sharex=True,figsize=(15,15))
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
    ylab_dict = {'Q': r'$Q_{in}\ [10^{3}\ m^{3}s^{-1}]$',
                'salt': r'$s_{in}\ [g\ kg^{-1}]$',
                'deltas': r'$\Delta s\ [g\ kg^{-1}]$'}
    # p_color = 'r'
    # m_color = 'b'
    lw = 2
        
    # fig = plt.figure()
    
    # ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
    # ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
    # ax3 = plt.subplot2grid((1,3), (0,2)) # map
    
    #ot = bulk['ot'] # (same as tef_df.index)
    
    axs[i].plot(ot,FR, color='tab:blue', linewidth=lw, label=r'$F_{R}$')
    axs[i].plot(ot,FE, color='tab:green', linewidth=lw, label=r'$F_{E}$')
    axs[i].plot(ot,FT, color='tab:red', linewidth=lw, label=r'$F_{T}$')
    axs[i].plot(ot,FTL, linestyle = '--', color='tab:pink', linewidth=lw, label=r'$F_{TL}$')
    axs[i].plot(ot,FTV, linestyle = ':', color='tab:orange', linewidth=lw, label=r'$F_{TV}$')
    axs[i].grid(True)    
    axs[i].set_ylabel(r'Salt flux $[m^{3}s^{-1} g\ kg^{-1}]$')
    axs[i].set_title(plot_title[i])
    axs[i].set_ylim(-4e4,4e4)
    
    
    
    # alpha=.3
    # ax1.plot(bulk['ot'],qp,'or',alpha=alpha)
    # ax1.plot(bulk['ot'],qm,'ob',alpha=alpha)
    
    # ax2.plot(bulk['ot'],sp,'or',alpha=alpha)
    # ax2.plot(bulk['ot'],sm,'ob',alpha=alpha)
    
    # ax2.plot(ot,tef_df['salt_p'].to_numpy()-tef_df['salt_m'].to_numpy(), color=plot_color[i], linewidth=lw, label=plot_label[i])
    # ax2.grid(True)
    # ax2.set_ylabel(ylab_dict['deltas'])
    # ax2.set_ylim(0,14)
    # ax2.set_xlim(pd.Timestamp('2020-04-01'), pd.Timestamp('2020-07-31'))

    # ax3.plot(ot,tef_df['salt_p'].to_numpy(), color=plot_color[i], linewidth=lw, label=plot_label[i])
    # ax3.grid(True)
    # ax3.set_ylabel(ylab_dict['salt'])
    # ax3.set_ylim(28,32)
    
    
    # # map
    # sn = sect_name.replace('.p','')
    # sinfo = sect_df.loc[sect_df.sn==sn,:]
    # i0 = sinfo.iloc[0,:].i
    # j0 = sinfo.iloc[0,:].j
    # uv0 = sinfo.iloc[0,:].uv
    # i1 = sinfo.iloc[-1,:].i
    # j1 = sinfo.iloc[-1,:].j
    # uv1 = sinfo.iloc[-1,:].uv
    # if uv0=='u':
    #     x0 = xu[j0,i0]
    #     y0 = yu[j0,i0]
    # elif uv0=='v':
    #     x0 = xv[j0,i0]
    #     y0 = yv[j0,i0]
    # if uv1=='u':
    #     x1 = xu[j1,i1]
    #     y1 = yu[j1,i1]
    # elif uv1=='v':
    #     x1 = xv[j1,i1]
    #     y1 = yv[j1,i1]
    # ax3.plot([x0,x1],[y0,y1],'-c', lw=3)
    # ax3.plot(x0,y0,'og', ms=10)
    # pfun.add_coast(ax3)
    # pfun.dar(ax3)
    # ax3.pcolormesh(xp, yp, -h, vmin=-100, vmax=100,
    #     cmap='jet', alpha=.4)
    
    # dx = x1-x0; dy = y1-y0
    # xmid = x0 + dx/2; ymid = y0 + dy/2
    # pad = np.max((np.sqrt(dx**2 + dy**2)*2,.1))
    # ax3.axis([x0-pad, x1+pad, y0-pad, y1+pad])
    # ax3.set_xlabel('Longitude [deg]')
    # ax3.set_ylabel('Latitude [deg]')
    # ax3.set_title(sn)
    # ax3.text(xmid - np.max((dy,.05))/6, ymid + np.max((dx,.05))/6, '+', fontweight='bold', c='r', fontsize=20,
    #     ha='center',va='center')
                
    # fig.tight_layout()
    
    # if Ldir['testing']:
    #     plt.show()
    # else:
    #     plt.savefig(out_dir / (sect_name.replace('.p','') + '.png'))
    #     plt.close()
axs[2].legend(loc='lower right')
axs[2].set_xlim(pd.Timestamp('2020-04-01'), pd.Timestamp('2020-07-31'))
plt.suptitle('Standard decomposition')
plt.savefig(out_dir / ('sd_plot_hourly.png'))
plt.close()
pfun.end_plot()
