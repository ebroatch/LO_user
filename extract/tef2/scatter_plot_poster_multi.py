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
in_dir = out_dir0 / ('bulk_' + Ldir['ds0'] + '_' + Ldir['ds1'])
out_dir = out_dir0 / ('bulk_poster_plots_' + Ldir['ds0'] + '_' + Ldir['ds1'])
Lfun.make_dir(out_dir, clean=True)

# sect_list = [item.name for item in in_dir.glob('*.p')]
# if Ldir['testing']:
#     sect_list = ['ss2.p']
sect_list_out = ['a1.p','a2.p','a3.p','a4.p','a5.p']
sect_list_sill = ['b1.p','b2.p','b3.p','b4.p','b5.p']
sect_list_in = ['c1.p','c2.p','c3.p','c4.p','c5.p']

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
plot_color_out = ['lightblue','tab:cyan','dodgerblue','tab:blue','blue']
plot_color_sill = ['gold','goldenrod','xkcd:yellow orange','tab:orange','peru']
plot_color_in = ['pink','tab:pink','mediumvioletred','tab:red','maroon']
plot_color = ['tab:purple','tab:blue','tab:green','tab:orange','tab:red']
#m_color = ['tab:cyan','xkcd:yellow orange','tab:pink']
plot_label_out = ['a1','a2','a3','a4','a5']
plot_label_sill = ['b1','b2','b3','b4','b5']
plot_label_in = ['c1','c2','c3','c4','c5']
fs = 12
plt.close('all')
pfun.start_plot(fs=fs, figsize=(21,10))

#fig, [ax1,ax2] = plt.subplots(2, 1, sharex=True,figsize=(15,10))
fig, [ax1,ax2,ax3] = plt.subplots(1, 3, figsize=(18,5))
# fig = plt.figure()   
# ax1 = plt.subplot2grid((2,3), (0,0), colspan=2) # Qin, Qout
# ax2 = plt.subplot2grid((2,3), (1,0), colspan=2) # Sin, Sout
# ax3 = plt.subplot2grid((1,3), (0,2)) # map

for i in range(len(sect_list_out)):
    sect_name = sect_list_out[i]
    bulk = pickle.load(open(in_dir / sect_name, 'rb'))

    tef_df = flux_fun.get_two_layer(in_dir, sect_name)
            
    # adjust units #leave for loglog plot 
    # tef_df['Q_p'] = tef_df['q_p']/1000
    # tef_df['Q_m'] = tef_df['q_m']/1000
    # tef_df['Q_prism'] = tef_df['qabs']/2000
    tef_df['Q_p'] = tef_df['q_p']
    tef_df['Q_m'] = tef_df['q_m']
    tef_df['Q_prism'] = tef_df['qabs']/2
                    
    lw = 2
    ot = bulk['ot'] # (same as tef_df.index)
    
    #ax1.scatter(tef_df['Q_prism'].to_numpy(),tef_df['Q_p'].to_numpy(), c=plot_color[i], linewidth=lw, label=plot_label[i])
    #ax1.plot(tef_df['Q_prism'].to_numpy(),tef_df['Q_p'].to_numpy(), color=plot_color[i], linewidth=lw, label=plot_label[i])
    ax1.loglog(tef_df['Q_prism'].to_numpy(),tef_df['Q_p'].to_numpy(), '-', lw=0.5, color=plot_color[i], label=plot_label_out[i])

for i in range(len(sect_list_sill)):
    sect_name = sect_list_sill[i]
    bulk = pickle.load(open(in_dir / sect_name, 'rb'))

    tef_df = flux_fun.get_two_layer(in_dir, sect_name)
            
    # adjust units #leave for loglog plot 
    # tef_df['Q_p'] = tef_df['q_p']/1000
    # tef_df['Q_m'] = tef_df['q_m']/1000
    # tef_df['Q_prism'] = tef_df['qabs']/2000
    tef_df['Q_p'] = tef_df['q_p']
    tef_df['Q_m'] = tef_df['q_m']
    tef_df['Q_prism'] = tef_df['qabs']/2
                    
    lw = 2
    ot = bulk['ot'] # (same as tef_df.index)
    
    #ax1.scatter(tef_df['Q_prism'].to_numpy(),tef_df['Q_p'].to_numpy(), c=plot_color[i], linewidth=lw, label=plot_label[i])
    #ax1.plot(tef_df['Q_prism'].to_numpy(),tef_df['Q_p'].to_numpy(), color=plot_color[i], linewidth=lw, label=plot_label[i])
    ax2.loglog(tef_df['Q_prism'].to_numpy(),tef_df['Q_p'].to_numpy(), '-', lw=0.5, color=plot_color[i], label=plot_label_sill[i])

for i in range(len(sect_list_in)):
    sect_name = sect_list_in[i]
    bulk = pickle.load(open(in_dir / sect_name, 'rb'))

    tef_df = flux_fun.get_two_layer(in_dir, sect_name)
            
    # adjust units #leave for loglog plot 
    # tef_df['Q_p'] = tef_df['q_p']/1000
    # tef_df['Q_m'] = tef_df['q_m']/1000
    # tef_df['Q_prism'] = tef_df['qabs']/2000
    tef_df['Q_p'] = tef_df['q_p']
    tef_df['Q_m'] = tef_df['q_m']
    tef_df['Q_prism'] = tef_df['qabs']/2
                    
    lw = 2
    ot = bulk['ot'] # (same as tef_df.index)
    
    #ax1.scatter(tef_df['Q_prism'].to_numpy(),tef_df['Q_p'].to_numpy(), c=plot_color[i], linewidth=lw, label=plot_label[i])
    #ax1.plot(tef_df['Q_prism'].to_numpy(),tef_df['Q_p'].to_numpy(), color=plot_color[i], linewidth=lw, label=plot_label[i])
    ax3.loglog(tef_df['Q_prism'].to_numpy(),tef_df['Q_p'].to_numpy(), '-', lw=0.5, color=plot_color[i], label=plot_label_in[i])

ax1.grid(True)
ax2.grid(True)
ax3.grid(True)
# ax1.axis('equal')
# ax2.axis('equal')
# ax3.axis('equal')
ax1.set_title('Outer basin sections')
ax2.set_title('Sill sections')
ax3.set_title('Inner basin sections')
# ax1.set_ylabel(r'$Q_{in} [10^{3}\ m^{3}s^{-1}]$')
# ax1.set_xlabel(r'$Q_{prism} [10^{3}\ m^{3}s^{-1}]$')
ax1.set_ylabel(r'$Q_{in} [m^{3}s^{-1}]$')
ax1.set_xlabel(r'$Q_{prism} [m^{3}s^{-1}]$')
ax2.set_ylabel(r'$Q_{in} [m^{3}s^{-1}]$')
ax2.set_xlabel(r'$Q_{prism} [m^{3}s^{-1}]$')
ax3.set_ylabel(r'$Q_{in} [m^{3}s^{-1}]$')
ax3.set_xlabel(r'$Q_{prism} [m^{3}s^{-1}]$')

ax1.set_xlim(left=1e4,right=5e4)
ax1.set_ylim(bottom=4e3, top=2e4)
ax2.set_xlim(left=1e4,right=3e4)
ax2.set_ylim(bottom=3e3, top=9e3)
ax3.set_xlim(left=8e2,right=4e4)
ax3.set_ylim(bottom=4e2, top=2e4)

# ax1.set_xlim(left=0)
ax1.legend(loc='lower right')
ax2.legend(loc='lower right')
ax3.legend(loc='lower right')
plt.savefig(out_dir / ('scatter_plot_poster_multi.png'))
plt.close()
pfun.end_plot()
