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
in_dir4 = out_dir0 / ('dvdk_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('sd_tef_plots_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

# sect_list = [item.name for item in in_dir.glob('*.p')]
# if Ldir['testing']:
#     sect_list = ['ss2.p']
# sect_list = ['a1.p','a3.p','b1.p','b3.p','b5.p','c3.p']
# sect_label = ['a1','a3','b1','b3','b5','c3']
sect_list = ['b1.p','b2.p','b3.p','b4.p','b5.p']
sect_nclist = ['b1.nc','b2.nc','b3.nc','b4.nc','b5.nc']
sect_label = ['b1','b2','b3','b4','b5']

# plot_color = ['tab:red','tab:orange','tab:green','tab:cyan','tab:blue','tab:purple']
plot_color = ['tab:cyan',plt.cm.Dark2(0),'tab:olive',plt.cm.Dark2(5),'tab:brown']
plot_color_light = [plt.cm.tab20(19),plt.cm.Set2(0),plt.cm.tab20(17),plt.cm.Set2(5),plt.cm.tab20(11)]

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

# PLOTTING
#plot_color = ['lightblue','tab:cyan','dodgerblue','tab:blue','blue','gold','goldenrod','xkcd:yellow orange','tab:orange','peru','pink','tab:pink','mediumvioletred','tab:red','maroon']
#plot_label = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
# p_color = ['tab:blue','tab:orange','tab:red']
# m_color = ['tab:cyan','xkcd:yellow orange','tab:pink']
# label_in = ['a3 in','b3 in','c3 in']
# label_out = ['a3 out','b3 out','c3 out']
# fs = 12
plt.close('all')
#pfun.start_plot(fs=fs, figsize=(21,10))
#pfun.start_plot(fs=fs, figsize=(7,10)) #narrower plot with one month xlim
pfun.start_plot(fs=14)


#fig, axs = plt.subplots(4, 1, sharex=True,figsize=(10,10),gridspec_kw={'height_ratios': [5,5,1,1]})
fig, axs = plt.subplots(2, 1, sharex=True,figsize=(8,8),gridspec_kw={'height_ratios': [1,8]})
# fig = plt.figure()   
# ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
# ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
# ax3 = plt.subplot2grid((1,3), (0,2)) # map

#Add Qprism and grey bars
axs[0].set_ylim(20,80)
axs[0].set_yticks([20,50,80])
# axs[1].set_ylim(-20,100) #to match when placed side by side
# axs[1].set_ylim(-5,65) #for FE, FE+FT, and tef
axs[1].set_ylim(18,42) #for only FE+FT and tef
# axs[2].set_ylim(-60,60)
# axs[3].set_ylim(-35,-25)

sect_name='b3'
bulk = xr.open_dataset(in_dir3 / (sect_name + '.nc'))
tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir3, sect_name)
tef_df['Q_p'] = tef_df['q_p']/1000
tef_df['Q_m'] = tef_df['q_m']/1000
tef_df['Q_prism']=tef_df['qprism']/1000
ot = bulk.time.values
lw=2
axs[0].plot(ot,tef_df['Q_prism'].to_numpy(), color='tab:gray', linewidth=lw)
snmid=(np.max(tef_df['Q_prism'].loc['2020-10-01':'2020-10-31'])+np.min(tef_df['Q_prism'].loc['2020-10-01':'2020-10-31']))/2
snbg=np.where(tef_df['Q_prism'].to_numpy()>snmid, 1, 0)
axs[0].pcolor(ot, axs[0].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
axs[1].pcolor(ot, axs[1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
# axs[2].pcolor(ot, axs[2].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
# axs[3].pcolor(ot, axs[3].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
axs[0].grid(True)

    # if sect_name=='b3.p':
    #     ds2 = xr.open_dataset(in_dir2 / sect_ncname) #changed to netcdf open in xarray
    #     qprism=zfun.lowpass(np.abs(ds2['qnet'].values-zfun.lowpass(ds2['qnet'].values, f='godin')), f='godin')/2 #for netcdf
    #     #qprism=zfun.lowpass(np.abs(ds2['qnet']-zfun.lowpass(ds2['qnet'], f='godin')), f='godin')/2
    #     pad=36
    #     qprism=(qprism[pad:-pad+1])[pad:-pad+1]

    #     axs[0].plot((ds2['time'].values[pad:-pad+1])[pad:-pad+1], qprism/1000, color='tab:grey', linewidth=lw)

for i in range(len(sect_list)):
    sect_name = sect_list[i]
    sect_ncname = sect_nclist[i]

    SD = xr.open_dataset(in_dir / sect_ncname)
    #SD = pickle.load(open(in_dir / sect_name, 'rb'))
    FR = SD['FR']
    FE = SD['FE']
    FT = SD['FT']
    FTL = SD['FTL']
    FTV = SD['FTV']
    F = SD['F']
    ot_sd = SD['time']

    # dvdk = xr.open_dataset(in_dir4 / sect_ncname)
    # F0 = dvdk['F0']
    # F1 = dvdk['F1']
    # F2 = dvdk['F2']
    # F2L = dvdk['F2L']
    # F2V = dvdk['F2V']
    # F_dvdk = dvdk['F']
    # ot_dvdk = dvdk['time']

    bulk = xr.open_dataset(in_dir3 / sect_ncname)
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir3, sect_label[i])
    tef_df['Q_p'] = tef_df['q_p']/1000
    tef_df['Q_m'] = tef_df['q_m']/1000
    tef_df['Q_prism']=tef_df['qprism']/1000
    ot_tef = bulk.time.values
                    
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
    # axs[3].plot(ot_dvdk, F0/1000, color=plot_color_light[i], linewidth=lw, label='$F_{0}$ '+sect_label[i])
    # axs[1].plot(ot_dvdk, F2/1000, color=plot_color_light[i], linewidth=lw, label='$F_{2}$ '+sect_label[i])
    # axs[2].plot(ot_dvdk, F1/1000, color=plot_color_light[i], linewidth=lw, label='$F_{1}$ '+sect_label[i])

    axs[1].plot(ot_tef,tef_df['Q_p'].to_numpy()*(tef_df['salt_p'].to_numpy()-tef_df['salt_m'].to_numpy()), color=plot_color_light[i], linewidth=lw, label='$Q_{in}\Delta s$ '+sect_label[i])

    # axs[3].plot(ot_sd, FR/1000, color=plot_color[i], linewidth=lw, ls='--', label='$F_{R}$ '+sect_label[i])
    axs[1].plot(ot_sd, (FE+FT)/1000, color=plot_color[i], linewidth=lw, ls='--', label='$F_{E}+F_{T}$ '+sect_label[i])
    # axs[1].plot(ot_sd, (FE)/1000, color=plot_color[i], linewidth=lw, ls=':', label='$F_{E}$ '+sect_label[i]) #OPTIONAL, COULD ADD BACK
    # axs[2].plot(ot_sd, FT/1000, color=plot_color[i], linewidth=lw, ls='--', label='$F_{T}$ '+sect_label[i])


    # axs[2].plot(ot, FTV, color=plot_color_light[i], linewidth=lw)#, ls=':')# label=sect_label[i])
    # axs[2].plot(ot, FTL, color=plot_color[i], linewidth=lw, ls=':')# label=sect_label[i])
    # axs[i].plot(ot,FTL, linestyle = '--', color='tab:pink', linewidth=lw, label=r'$F_{TL}$')
    # axs[i].plot(ot,FTV, linestyle = ':', color='tab:orange', linewidth=lw, label=r'$F_{TV}$')

axs[0].grid(True)
axs[1].grid(True) 
# axs[2].grid(True)
# axs[3].grid(True)

axs[0].text(0.02,0.8,'A',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[0].transAxes)
axs[1].text(0.02,0.95,'B',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[1].transAxes)
# axs[2].text(0.02,0.95,'C',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[2].transAxes)
# axs[3].text(0.02,0.8,'D',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[3].transAxes)

axs[1].legend(loc='upper right',fontsize=12)
# axs[2].legend(loc='upper right',fontsize=12)

# axs[2].set_xlim(pd.Timestamp('2020-10-01'), pd.Timestamp('2020-11-15'))
# axs[3].set_xticks(ticks=[pd.Timestamp('2020-10-01'), pd.Timestamp('2020-10-15'), pd.Timestamp('2020-11-01'), pd.Timestamp('2020-11-15')])
#axs[3].set_xlim(pd.Timestamp('2020-09-01'), pd.Timestamp('2020-12-31'))
axs[1].set_xlim(pd.Timestamp('2020-10-01'), pd.Timestamp('2020-10-31'))
#plt.xticks(rotation=90)

# axs[3].set_ylabel('$F_{0}$ or $F_{R}$\n$[10^{3}\ m^{3}s^{-1} g\ kg^{-1}]$')
axs[1].set_ylabel('Salt transport\n$[10^{3}\ m^{3}s^{-1} g\ kg^{-1}]$')
# axs[2].set_ylabel('$F_{1}$ or $F_{T}$\n$[10^{3}\ m^{3}s^{-1} g\ kg^{-1}]$')
# axs[0].set_ylabel('$Q_{prism}$\n$[10^{3}\ m^{3}s^{-1}]$')
axs[0].set_ylabel('$Q_{prism}$ (b3)\n$[10^{3}\ m^{3}s^{-1}]$')

axs[1].xaxis.set_major_formatter(mdates.DateFormatter('%j')) #yearday
axs[1].set_xlabel('Yearday')

# if Ldir['gridname']=='sill20kmdeep':
#     axs[0].set_title('20km sill')
# elif Ldir['gridname']=='sill5km':
#     axs[0].set_title('5km sill')
# elif Ldir['gridname']=='sill80km':
#     axs[0].set_title('80km sill')
#plt.suptitle('Standard decomposition')
# axs[0].set_title(Ldir['gtagex'])
plt.tight_layout()
plt.savefig(out_dir / ('sd_tef_plot_hourly.png'))
plt.close()
pfun.end_plot()
    
    
    
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

# #scatter plot
# fig, ax = plt.subplots(1, 1, figsize=(10,10))
# for i in range(len(sect_list)):
#     sect_name = sect_list[i]
#     sect_ncname = sect_nclist[i] #change to netcdf
    
#     SD = xr.open_dataset(in_dir / sect_ncname)
#     # FR = SD['FR']
#     FE = SD['FE']
#     FT = SD['FT']
#     # FTL = SD['FTL']
#     # FTV = SD['FTV']
#     # F = SD['F']
#     # ot = SD['time']

#     #ds2= pickle.load(open(in_dir2 / sect_name, 'rb'))
#     #qprism=zfun.lowpass(np.abs(ds2['qnet']-zfun.lowpass(ds2['qnet'], f='godin')), f='godin')/2
#     ds2 = xr.open_dataset(in_dir2 / sect_ncname) #changed to netcdf open in xarray
#     qprism=zfun.lowpass(np.abs(ds2['qnet'].values-zfun.lowpass(ds2['qnet'].values, f='godin')), f='godin')/2 #for netcdf
#     pad=36
#     qprism=(qprism[pad:-pad+1])[pad:-pad+1]

#     ax.plot(qprism, FT/(FE+FT), color=plot_color[i], linewidth=lw, label=sect_label[i])

# ax.grid(True)
# ax.set_ylabel('FT/(FE+FT)')
# ax.set_xlabel('Qprism')
# ax.legend(loc='lower right')
# ax.set_title(Ldir['gtagex'])

# plt.savefig(out_dir / ('sd_plot_scatter.png'))
# plt.close()
