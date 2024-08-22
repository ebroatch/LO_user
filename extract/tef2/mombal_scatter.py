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

gridnamelist = ['sill5km','sill10km','sill20kmdeep','sill40km','sill80km']
taglist = ['t0','t2','t2','t2','t2']
gtaglist = ['sill5km_t0','sill10km_t2','sill20kmdeep_t2','sill40km_t2','sill80km_t2']
gtagexlist = ['sill5km_t0_xa0','sill10km_t2_xa0','sill20kmdeep_t2_xa0','sill40km_t2_xa0','sill80km_t2_xa0']
gridlist = [Ldir['data'] / 'grids' / 'sill5km', Ldir['data'] / 'grids' / 'sill10km', Ldir['data'] / 'grids' / 'sill20kmdeep', Ldir['data'] / 'grids' / 'sill40km', Ldir['data'] / 'grids' / 'sill80km']

silllenstr = ['5km','10km','20km','40km','80km']
silllen = [5e3,10e3,20e3,40e3,80e3]

fig, ax = plt.subplots(1, 1, figsize=(8,8))#,gridspec_kw={'height_ratios': [3,1,1,1]})#figsize 20,10 for 3 sects

for gi in range(len(gtagexlist)):
    Ldir['gridname']=gridnamelist[gi]
    Ldir['tag']=taglist[gi]
    Ldir['gtag']=gtaglist[gi]
    Ldir['gtagex']=gtagexlist[gi]
    Ldir['grid']=gridlist[gi]

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

    #sect_list = ['b1','b2','b3','b4','b5']
    sect_list = ['b1','b3','b5'] #only 3 sections for readability
    #plot_label = ['b1','b2','b3','b4','b5']
    #plot_color = ['tab:red','tab:orange','tab:green','tab:cyan','tab:blue']
    plot_marker = ['^','o','s']

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
    # plt.close('all')
    #pfun.start_plot(fs=fs, figsize=(10,20))

    #fig, [ax1,ax2,ax3] = plt.subplots(3, 1, sharex=True,figsize=(15,15))
    # fig, [ax0,ax1,ax2,ax3] = plt.subplots(4, 1, sharex=True,figsize=(15,7.7),gridspec_kw={'height_ratios': [1,4,2,2]})
    #fig, axs = plt.subplots(len(sect_list), 2, sharex=True, figsize=(20,15))#,gridspec_kw={'height_ratios': [1,4,2,2,2,2]})#figsize 20,10 for 3 sects
    
    # fig, axs = plt.subplots(4, 1, figsize=(8,12),gridspec_kw={'height_ratios': [3,1,1,1]})#figsize 20,10 for 3 sects
    
    # fig = plt.figure()   
    # ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
    # ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
    # ax3 = plt.subplot2grid((1,3), (0,2)) # map

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

        # get tide info from the tide excursion calculator
        excur_dir = out_dir0 / ('tide_excursion_' + Ldir['ds0'] + '_' + Ldir['ds1'])
        te_fn = excur_dir / ('TE_b3.p') #could change if using other sections
        TE = pd.read_pickle(te_fn)

                        
        # labels and colors
        # ylab_dict = {'Q': r'Transport $[10^{3}\ m^{3}s^{-1}]$',
        #             'salt': r'Salinity $[g\ kg^{-1}]$'}
        ylab_dict = {'Qprism': '$Q_{prism}$\n$[10^{3}\ m^{3}s^{-1}]$',
                    'Q': '$Q_{in}$\n$[10^{3}\ m^{3}s^{-1}]$',
                    'sin': '$s_{in}$\n$[g\ kg^{-1}]$',
                    'sout': '$s_{out}$\n$[g\ kg^{-1}]$',
                    'deltas': '$\Delta s$\n$[g\ kg^{-1}]$'}
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
        # ax.scatter(i,dudt_in_alt,color='tab:red', label='d/dt(Qin/Ain)',marker=plot_marker[i])
        # axs[i,0].plot(i,tef_df['coriolis_p'],color='tab:purple', label='coriolis')
        ax.errorbar(silllen[i]/TE['TE_spring'],np.mean(tef_df['pg_p']),yerr=np.std(tef_df['pg_p']),color='tab:green', label='PG',marker=plot_marker[i],markersize=10,lw=0.5,capsize=6,ls=None)
        ax.errorbar(silllen[i]/TE['TE_spring'],np.mean(tef_df['stressdiv_p']),yerr=np.std(tef_df['stressdiv_p']),color='tab:blue', label='stressdiv',marker=plot_marker[i],markersize=10,lw=0.5,capsize=6,ls=None)
        resid=tef_df['dudt_p']-tef_df['coriolis_p']-tef_df['pg_p']-tef_df['stressdiv_p']
        ax.errorbar(silllen[i]/TE['TE_spring'],np.mean(resid),yerr=np.std(resid),color='k', label='residual (advection)',marker=plot_marker[i],markersize=10,lw=0.5,capsize=6,ls=None)
        if gi==0:
            if i==0:
                ax.legend(loc='lower right')

        # axs[i,1].plot(ot,tef_df['dudt_m'],color='tab:red',ls='--', label='du/dt')
        # axs[i,1].plot(ot,dudt_out_alt,color='tab:red',ls='--', label='d/dt(Qout/Aout)')
        # axs[i,1].plot(ot,tef_df['coriolis_m'],color='tab:purple', ls='--', label='coriolis')
        # axs[i,1].plot(ot,tef_df['pg_m'],color='tab:green', ls='--', label='PG')
        # axs[i,1].plot(ot,tef_df['stressdiv_m'],color='tab:blue', ls='--', label='stressdiv')
        # axs[i,1].plot(ot,tef_df['dudt_m']-tef_df['coriolis_m']-tef_df['pg_m']-tef_df['stressdiv_m'],ls='--',color='k', label='residual (advection)')
        #axs[i].plot(ot,tef_df['dudt_p']-tef_df['coriolis_p']-tef_df['pg_p']-tef_df['stressdiv_p'],color='k', label='residual (advection) in')
        # axs[i,0].text(0.05,0.9,sect_name+' in',transform=axs[i,0].transAxes)
        # axs[i,1].text(0.05,0.9,sect_name+' out',transform=axs[i,1].transAxes)
        
        # axs[1].grid(True)
        # axs[2].grid(True)
        # axs[3].grid(True)
        # if i==1:
        #     axs[1].plot(ot,Qin,color='tab:red', label='Qin')
        #     axs[1].plot(ot,Qout,color='tab:blue', label='Qin')
        #     axs[2].plot(ot,Uin,color='tab:red', label='Uin')
        #     axs[2].plot(ot,Uout,color='tab:blue', label='Uout')
        #     axs[3].plot(ot,Ain,color='tab:red', label='Ain')
        #     axs[3].plot(ot,Aout,color='tab:blue', label='Aout')
        #     axs[1].set_title('20km b3 values')
        #     axs[1].legend()
        #     axs[2].legend()
        #     axs[3].legend()

        
        #can reset these axes later
        # axs[i,0].set_ylim(-0.0005,0.0005)
        # axs[i,1].set_ylim(-0.0005,0.0005)
        # if i==0:
        #     axs[i,0].set_ylim(-0.004,0.004)
        #     # axs[i,1].set_ylim(-0.001,0.001)
        # if i==4:
        #     axs[i,0].set_ylim(-0.001,0.001)
        #     axs[i,1].set_ylim(-0.002,0.002)


ax.grid(True)
ax.set_xlabel('Lsill/Ltide (spring)')
ax.set_ylabel('Momentum balance terms')
# axs[1].set_xlabel('Yearday')
ax.text(0.05,0.9,'triangle=b1,\ncircle=b3\nsquare=b5',transform=ax.transAxes)
#axs[1].legend(loc='lower right')
#ax.legend(loc='lower right')
#plt.suptitle(Ldir['gtagex']+'outflow hourly flux-weighted momentum balance')
plt.suptitle('Momentum balance for inflow (mean&stdev)')
plt.savefig(out_dir / ('mombal_scatter.png'))
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