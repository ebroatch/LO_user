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
import tef_fun

from lo_tools import extract_argfun as exfun
Ldir = exfun.intro() # this handles the argument passing

gctag = Ldir['gridname'] + '_' + Ldir['collection_tag']
tef2_dir = Ldir['LOo'] / 'extract' / 'tef2'

sect_df_fn = tef2_dir / ('sect_df_' + gctag + '.p')
sect_df = pd.read_pickle(sect_df_fn)

out_dir0 = Ldir['LOo'] / 'extract' / Ldir['gtagex'] / 'tef2'
in_dir = out_dir0 / ('bulk_mombal_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
in_dir2 = out_dir0 / ('bulk_mombal_area_hourly_' + Ldir['ds0'] + '_' + Ldir['ds1'])
# in_dir = out_dir0 / ('bulk_mombal_' + Ldir['ds0'] + '_' + Ldir['ds1'])
# in_dir2 = out_dir0 / ('bulk_mombal_area_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('bulk_mombal_plots_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

# sect_list = [item.name for item in in_dir.glob('*.p')]
# if Ldir['testing']:
#     sect_list = ['ss2.p']
#sect_list = ['a3.p','b3.p','c3.p']
#sect_list = ['a1.p','a2.p','a3.p','a4.p','a5.p','b1.p','b2.p','b3.p','b4.p','b5.p','c1.p','c2.p','c3.p','c4.p','c5.p']
#sect_list = ['a1.p','a3.p','b1.p','b2.p','b3.p','b4.p','b5.p','c3.p']
#plot_label = ['a1','a3','b1','b2','b3','b4','b5','c3']
#plot_color = ['tab:red','tab:orange','tab:olive','tab:green','tab:cyan','tab:blue','tab:purple','tab:pink']
#plot_color = ['k','tab:gray','tab:red','tab:orange','tab:green','tab:cyan','tab:blue','tab:brown']

# sect_list = ['b1','b2','b3','b4','b5']
sect_list = ['b1','b3','b5'] #only 3 sections for readability
#plot_label = ['b1','b2','b3','b4','b5']
# plot_color = ['tab:red','tab:orange','tab:green','tab:cyan','tab:blue']

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
#pfun.start_plot(fs=fs, figsize=(10,20))

pfun.start_plot(fs=14)
#fig, [ax1,ax2,ax3] = plt.subplots(3, 1, sharex=True,figsize=(15,15))
# fig, [ax0,ax1,ax2,ax3] = plt.subplots(4, 1, sharex=True,figsize=(15,7.7),gridspec_kw={'height_ratios': [1,4,2,2]})
# fig, axs = plt.subplots(len(sect_list)+1, 1, sharex=True, figsize=(8,12))#,gridspec_kw={'height_ratios': [1,4,2,2,2,2]})#figsize 20,10 for 3 sects
fig, axs = plt.subplots(len(sect_list)+1, 1, sharex=True, figsize=(8,10),gridspec_kw={'height_ratios': [1,4,4,4]})#figsize 20,10 for 3 sects
# fig = plt.figure()   
# ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
# ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
# ax3 = plt.subplot2grid((1,3), (0,2)) # map

#Add Qprism and grey bars
sect_name='b3'
bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
tef_df['Q_p'] = tef_df['q_p']/1000
tef_df['Q_m'] = tef_df['q_m']/1000
tef_df['Q_prism']=tef_df['qprism']/1000
ot = bulk.time.values
lw=2
axs[0].plot(ot,tef_df['Q_prism'].to_numpy(), color='tab:gray', linewidth=lw)
axs[0].set_ylabel('$Q_{prism}$ (b3)\n$[10^{3}\ m^{3}s^{-1}]$')
axs[0].set_ylim(20,100)
snmid=(np.max(tef_df['Q_prism'].loc['2020-10-01':'2020-10-31'])+np.min(tef_df['Q_prism'].loc['2020-10-01':'2020-10-31']))/2
snbg=np.where(tef_df['Q_prism'].to_numpy()>snmid, 1, 0)
axs[0].pcolor(ot, axs[0].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
axs[1].pcolor(ot, axs[1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
axs[2].pcolor(ot, axs[2].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
axs[3].pcolor(ot, axs[3].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
axs[0].grid(True)

#Plot mombal terms
for i in range(len(sect_list)):
    sect_name = sect_list[i]
    bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
            
    # adjust units
    tef_df['Q_p'] = tef_df['q_p']/1000
    tef_df['Q_m'] = tef_df['q_m']/1000
    tef_df['Q_prism']=tef_df['qprism']/1000

    # bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))
    tef_df_area, vn_list_area, vec_list_area = tef_fun.get_two_layer(in_dir2, sect_name)
    Qin = tef_df_area['q_p']
    Qout = tef_df_area['q_m']
    Ain = tef_df_area['a_p']
    Aout = tef_df_area['a_m']
    Uin = Qin/Ain
    Uout = Qout/Aout
    dudt_in_alt =np.concatenate(([np.nan],(Uin.values[2:]-Uin.values[:-2])/(2*3600),[np.nan])) #REMOVE 24 FOR HOURLY DATA
    dudt_out_alt =np.concatenate(([np.nan],(Uout.values[2:]-Uout.values[:-2])/(2*3600),[np.nan]))
              
    # labels and colors
    # ylab_dict = {'Q': r'Transport $[10^{3}\ m^{3}s^{-1}]$',
    #             'salt': r'Salinity $[g\ kg^{-1}]$'}
    # ylab_dict = {'Qprism': '$Q_{prism}$\n$[10^{3}\ m^{3}s^{-1}]$',
    #             'Q': '$Q_{in}$\n$[10^{3}\ m^{3}s^{-1}]$',
    #             'sin': '$s_{in}$\n$[g\ kg^{-1}]$',
    #             'sout': '$s_{out}$\n$[g\ kg^{-1}]$',
    #             'deltas': '$\Delta s$\n$[g\ kg^{-1}]$'}
    # p_color = 'r'
    # m_color = 'b'
    lw = 2
        
    # fig = plt.figure()
    
    # ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
    # ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
    # ax3 = plt.subplot2grid((1,3), (0,2)) # map

    #ot = bulk['ot'] # (same as tef_df.index)
    ot = bulk.time.values
    
    # axs[i,0].plot(ot,tef_df['dudt_p'],color='tab:red', label='du/dt')
    # axs[i].plot(ot,dudt_in_alt,color='tab:olive', label='d/dt(Qin/Ain)')
    axs[i+1].plot(ot,dudt_in_alt,color='tab:olive', label='Acceleration')
    axs[i+1].plot(ot,tef_df['coriolis_p'],color='tab:purple', label='Coriolis')
    axs[i+1].plot(ot,tef_df['pg_p'],color='tab:red', label='Pressure gradient')
    axs[i+1].plot(ot,tef_df['stressdiv_p'],color='tab:cyan', label='Stress divergence')
    axs[i+1].plot(ot,tef_df['dudt_p']-tef_df['coriolis_p']-tef_df['pg_p']-tef_df['stressdiv_p'],color='k', label='Residual (advection)')
    if i==2:
        axs[i+1].legend(loc='upper center',fontsize=12,bbox_to_anchor=(0.5, -0.05), ncol=5, mode="expand", borderaxespad=0.)
    # axs[i,1].plot(ot,tef_df['dudt_m'],color='tab:red',ls='--', label='du/dt')
    axs[i+1].plot(ot,dudt_out_alt,color='tab:olive',ls='--', label='d/dt(Qout/Aout)')
    axs[i+1].plot(ot,tef_df['coriolis_m'],color='tab:purple', ls='--', label='coriolis')
    axs[i+1].plot(ot,tef_df['pg_m'],color='tab:red', ls='--', label='PG')
    axs[i+1].plot(ot,tef_df['stressdiv_m'],color='tab:cyan', ls='--', label='stressdiv')
    axs[i+1].plot(ot,tef_df['dudt_m']-tef_df['coriolis_m']-tef_df['pg_m']-tef_df['stressdiv_m'],ls='--',color='k', label='residual (advection)')
    #axs[i].plot(ot,tef_df['dudt_p']-tef_df['coriolis_p']-tef_df['pg_p']-tef_df['stressdiv_p'],color='k', label='residual (advection) in')
    # axs[i].text(0.05,0.9,sect_name,transform=axs[i].transAxes)
    # axs[i,0].text(0.05,0.9,sect_name+' in',transform=axs[i,0].transAxes)
    # axs[i,1].text(0.05,0.9,sect_name+' out',transform=axs[i,1].transAxes)
    axs[i+1].grid(True)
    # axs[i,1].grid(True)
    
    if Ldir['gridname']=='sill5km':
        ylimmax=0.0025
    else:
        ylimmax=0.004
    axs[i+1].set_ylim(-ylimmax,ylimmax)

    # axs[i,0].set_ylim(-0.0005,0.0005)
    # axs[i,1].set_ylim(-0.0005,0.0005)
    # if i==0:
    #     axs[i,0].set_ylim(-0.004,0.004)
    #     # axs[i,1].set_ylim(-0.001,0.001)
    # if i==4:
    #     axs[i,0].set_ylim(-0.001,0.001)
    #     axs[i,1].set_ylim(-0.002,0.002)


    # ax1.plot(ot,tef_df['Q_p'].to_numpy(), color=plot_color[i], linewidth=lw, label=sect_name)
    # #ax1.plot(ot,tef_df['Q_m'].to_numpy(), color=m_color[i], linewidth=lw, label=label_out[i])
    # ax1.grid(True)
    # ax1.set_ylabel(ylab_dict['Q'])
    # ax1.set_ylim(0,20)
    #ax1.set_ylim(0,16)
    #ax1.set_yticks(ticks=[0,4,8,12,16])
    
    # qp = bulk['q'].copy()/1000
    # qp[qp<0] = np.nan
    # qm = bulk['q'].copy()/1000
    # qm[qm>0]=np.nan
    # sp = bulk['salt'].copy()
    # sp[np.isnan(qp)] = np.nan
    # sm = bulk['salt'].copy()
    # sm[np.isnan(qm)]=np.nan
    
    # alpha=.3
    # ax1.plot(bulk['ot'],qp,'or',alpha=alpha)
    # ax1.plot(bulk['ot'],qm,'ob',alpha=alpha)
    
    # ax2.plot(bulk['ot'],sp,'or',alpha=alpha)
    # ax2.plot(bulk['ot'],sm,'ob',alpha=alpha)

    # ax2.plot(ot,tef_df['salt_p'].to_numpy(), color=plot_color[i], linewidth=lw, label=sect_name)
    # ax2.grid(True)
    # ax2.set_ylabel(ylab_dict['sin'])
    # ax2.set_ylim(24,34)

    # ax3.plot(ot,tef_df['salt_m'].to_numpy(), color=plot_color[i], linewidth=lw, label=sect_name)
    # ax3.grid(True)
    # ax3.set_ylabel(ylab_dict['sout'])
    # ax3.set_ylim(22,32)
    
    # ax4.plot(ot,tef_df['salt_p'].to_numpy()-tef_df['salt_m'].to_numpy(), color=plot_color[i], linewidth=lw, label=sect_name)
    # ax4.grid(True)
    # ax4.set_ylabel(ylab_dict['deltas'])
    # ax4.set_ylim(0,10)
    # #ax2.set_ylim(0,10)
    # #ax2.set_xlim(pd.Timestamp('2020-09-01'), pd.Timestamp('2020-12-31'))

    # ax5.plot(ot,tef_df['Q_p'].to_numpy()*(tef_df['salt_p'].to_numpy()-tef_df['salt_m'].to_numpy()), color=plot_color[i], linewidth=lw, label=sect_name)
    # ax5.grid(True)
    # ax5.set_ylabel('$Q_{in}\Delta s$\n$[10^{3}\ m^{3}s^{-1} g\ kg^{-1}]$')
    # ax5.set_ylim(0,75)

    #ax4.set_xlim(pd.Timestamp('2020-09-01'), pd.Timestamp('2020-12-31'))


    # if i==0:
    #     ax0.plot(ot,tef_df['Q_prism'].to_numpy(), color='tab:gray', linewidth=lw)
    #     ax0.set_ylabel(ylab_dict['Qprism'])
    #     ax0.set_ylim(20,100)
    #     #ax0.set_yticks(ticks=[20,30,40,50])
    #     # ax0.set_xlim(pd.Timestamp('2020-09-01'), pd.Timestamp('2020-12-31'))
    #     snmid=(np.max(tef_df['Q_prism'].loc['2020-10-01':'2020-10-31'])+np.min(tef_df['Q_prism'].loc['2020-10-01':'2020-10-31']))/2
    #     snbg=np.where(tef_df['Q_prism'].to_numpy()>snmid, 1, 0)
    #     ax0.pcolor(ot, ax0.get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
    #     ax1.pcolor(ot, ax1.get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
    #     ax2.pcolor(ot, ax2.get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
    #     ax3.pcolor(ot, ax3.get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
    #     ax4.pcolor(ot, ax4.get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
    #     ax5.pcolor(ot, ax5.get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=0, vmax=2, alpha=0.3, linewidth=0, antialiased=True)
    #     ax0.grid(True)

    
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

# axs[0].text(0.05,0.9,'b1',fontweight='bold',fontsize=14,transform=axs[0].transAxes)
# axs[1].text(0.05,0.9,'b3',fontweight='bold',fontsize=14,transform=axs[1].transAxes)
# axs[2].text(0.05,0.9,'b5',fontweight='bold',fontsize=14,transform=axs[2].transAxes)

axs[0].text(0.02,0.8,'A',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[0].transAxes)
axs[1].text(0.02,0.95,'B',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[1].transAxes)
axs[2].text(0.02,0.95,'C',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[2].transAxes)
axs[3].text(0.02,0.95,'D',ha='left',va='top',fontweight='bold',fontsize=18,transform=axs[3].transAxes)
axs[1].text(0.99,0.98,'section b1',fontsize=10,ha='right',va='top',transform=axs[1].transAxes)
axs[2].text(0.99,0.98,'section b3',fontsize=10,ha='right',va='top',transform=axs[2].transAxes)
axs[3].text(0.99,0.98,'section b5',fontsize=10,ha='right',va='top',transform=axs[3].transAxes)

axs[3].set_xlim(pd.Timestamp('2020-10-01'), pd.Timestamp('2020-10-31'))
# axs[4,1].set_xlim(pd.Timestamp('2020-10-01'), pd.Timestamp('2020-11-30'))
# axs[2,0].xaxis.set_major_formatter(mdates.DateFormatter('%-d'))
# axs[2,1].xaxis.set_major_formatter(mdates.DateFormatter('%-d'))
axs[3].xaxis.set_major_formatter(mdates.DateFormatter('%j')) #yearday
# axs[4,1].xaxis.set_major_formatter(mdates.DateFormatter('%j'))
# axs[2,0].set_xlabel('Day')
# axs[2,1].set_xlabel('Day')
axs[1].set_ylabel('Momentum balance terms\n$[m\ s^{-2}]$')
axs[2].set_ylabel('Momentum balance terms\n$[m\ s^{-2}]$')
axs[3].set_ylabel('Momentum balance terms\n$[m\ s^{-2}]$')
axs[3].set_xlabel('Yearday')
# axs[4,1].set_xlabel('Yearday')
# axs[4,0].legend(loc='lower right')
# axs[4,1].legend(loc='lower right')
#axs[0].legend(loc='lower right')
# plt.suptitle(Ldir['gtagex']+' hourly flux-weighted momentum balance')
plt.tight_layout()
plt.savefig(out_dir / ('bulk_mombal_plot3.png'))
plt.close()
pfun.end_plot()

# #Scatter plot
# fig, [ax1,ax2,ax3] = plt.subplots(1, 3 ,figsize=(18,7))
# for i in range(len(sect_list)):
#     sect_name = sect_list[i]
#     bulk = xr.open_dataset(in_dir / (sect_name + '.nc'))

#     tef_df, vn_list, vec_list = tef_fun.get_two_layer(in_dir, sect_name)
            
#     # adjust units
#     tef_df['Q_p'] = tef_df['q_p']/1000
#     tef_df['Q_m'] = tef_df['q_m']/1000
#     tef_df['Q_prism']=tef_df['qprism']/1000
                    
#     # labels and colors
#     # ylab_dict = {'Q': r'Transport $[10^{3}\ m^{3}s^{-1}]$',
#     #             'salt': r'Salinity $[g\ kg^{-1}]$'}
#     ylab_dict = {'Qprism': '$Q_{prism}$\n$[10^{3}\ m^{3}s^{-1}]$',
#                 'Q': '$Q_{in}$\n$[10^{3}\ m^{3}s^{-1}]$',
#                 'salt': '$s_{in}$\n$[g\ kg^{-1}]$',
#                 'deltas': '$\Delta s$\n$[g\ kg^{-1}]$'}

#     pad=36
#     ax1.plot(tef_df['Q_prism'][pad:-pad+1].to_numpy(),tef_df['Q_p'][pad:-pad+1].to_numpy(), '-', lw=0.5, color=plot_color[i], label=sect_name) #cut out first couple of days for weird qprism
#     ax2.plot(tef_df['Q_prism'][pad:-pad+1].to_numpy(),tef_df['salt_p'][pad:-pad+1].to_numpy()-tef_df['salt_m'][pad:-pad+1].to_numpy(), '-', lw=0.5, color=plot_color[i], label=sect_name)
#     ax3.plot(tef_df['Q_prism'][pad:-pad+1].to_numpy(),tef_df['Q_p'][pad:-pad+1].to_numpy()*(tef_df['salt_p'][pad:-pad+1].to_numpy()-tef_df['salt_m'][pad:-pad+1].to_numpy()), '-', lw=0.5, color=plot_color[i], label=sect_name)

# ax1.set_ylabel(r'$Q_{in} [m^{3}s^{-1}]$')
# ax1.set_xlabel(r'$Q_{prism} [m^{3}s^{-1}]$')
# ax2.set_ylabel(r'$\Delta s [g kg^{-1}]$')
# ax2.set_xlabel(r'$Q_{prism} [m^{3}s^{-1}]$')
# ax3.set_ylabel(r'$Q_{in} \Delta s [m^{3}s^{-1}g kg^{-1}]$')
# ax3.set_xlabel(r'$Q_{prism} [m^{3}s^{-1}]$')

# ax1.set_box_aspect(1)
# ax2.set_box_aspect(1)
# ax3.set_box_aspect(1)

# ax1.grid(True)
# ax2.grid(True)
# ax3.grid(True)

# ax3.legend(loc='lower right')
# ax2.set_title(Ldir['gtagex'])
# plt.savefig(out_dir / ('tef_plot_scatter.png'))
# plt.close()
# pfun.end_plot()