"""
Module of plotting functions.

Each function creates, and optionally saves, a plot of fields
from a ROMS history file.

INPUT: in_dict: a tuple with information to pass to the plot, such as:
- fn: text string with the full path name of the history file to plot
- fn_out: text string with full path of output file name
- auto_vlims: a boolean governing how color limits are set
- testing: a boolean for testing (e.g. shorter, faster particle tracking)
OUTPUT: either a screen image or a graphics file

"""
import numpy as np
import xarray as xr
import pickle
from datetime import datetime, timedelta
import pandas as pd
from cmocean import cm
import copy

from lo_tools import Lfun, zfun, zrfun
from lo_tools import plotting_functions as pfun
import pinfo
from importlib import reload
reload(pfun)
reload(pinfo)

Ldir = Lfun.Lstart()
if '_mac' in Ldir['lo_env']: # mac version
    pass
else: # remote linux version
    import matplotlib as mpl
    mpl.use('Agg')
import matplotlib.pyplot as plt
    

def P_sect_eb(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'salt'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    lon = G['lon_rho']
    lat = G['lat_rho']
    zdeep = -205
    #x = np.linspace(1.1, -1, 500) #sill1
    x = np.linspace(1.2, -1, 500) #sill2
    xcoast = 94.3 #sill2
    y = 45 * np.ones(x.shape)
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING
    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    # fig.colorbar(cs, ax=ax) # It is identical to that of the section
    #pfun.add_coast(ax)
    aaf = [-4, 2, 43, 47] # focus domain
    ax.axis(aaf)
    pfun.dar(ax)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'r', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    ax.set_xticks([-4, -2, 0, 2])
    ax.set_yticks([43, 45, 47])
    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], 'k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-k', linewidth=2)
    ax.set_xlim(dist.min(), dist.max())
    ax.invert_xaxis()
    ax.set_xticks([0, 40, 80, 120, 160, 200])
    ax.vlines(xcoast, -200, 0, linestyles='dashed') #add coast line
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:],
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_zoom_eb(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    Uses the new pfun.get_sect() function.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])

    # PLOT CODE
    vn = 'salt'
    if vn == 'salt':
        pinfo.cmap_dict[vn] = 'Spectral_r'

    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])

    # CREATE THE SECTION
    # create track by hand
    lon = G['lon_rho']
    lat = G['lat_rho']
    zdeep = -205
    #x = np.linspace(1.1, -1, 500) #sill1
    x_e = np.linspace(0.45, 0.675, 100) #use less points for shorter section?
    #xcoast = 94.3 #sill2
    y_e = 45 * np.ones(x_e.shape)

    #v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict) #old section function
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = pfun.get_sect(in_dict['fn'], vn, x_e, y_e)

    # COLOR
    # scaled section data
    #sf = pinfo.fac_dict[vn] * v3['sectvarf'] #old
    sf = pinfo.fac_dict[vn] * fld_s #new

    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING
    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    # fig.colorbar(cs, ax=ax) # It is identical to that of the section
    #pfun.add_coast(ax)
    aaf = [0.45, 0.675, 44.95, 45.05] # focus domain
    #aaf = [-4, 2, 43, 47] # focus domain
    ax.axis(aaf)
    pfun.dar(ax)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    # ax.plot(x[idist0], y[idist0], 'r', markersize=5, markerfacecolor='w',
    #     markeredgecolor='r', markeredgewidth=2) #old
    ax.plot(x[0], y[0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2) #new
    ax.set_xticks([0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675])
    ax.set_yticks([44.95, 45, 45.05])
    # ax.set_xticks([-4, -2, 0, 2])
    # ax.set_yticks([43, 45, 47])

    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    # ax.plot(dist, v2['zbot'], 'k', linewidth=2)
    # ax.plot(dist, v2['zeta'], '-k', linewidth=2)
    ax.plot(dist, zbot, '-k', linewidth=2)
    ax.plot(dist, ztop, '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    #ax.invert_xaxis() #not sure if this needs to be changed?
    #ax.set_xticks([0, 2, 4, 6, 8, 10, 12]) #comment out for now
    #ax.vlines(xcoast, -200, 0, linestyles='dashed') #add coast line #don't need since zoomed in
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    # cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:],
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #old
    cs = ax.pcolormesh(dist_se,zw_se,sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #new
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_contour_eb(in_dict):
    """
    This plots a section (distance, z) with line contours, and makes sure
    that the color limits are identical.
    
    Uses the new pfun.get_sect() function.
    """
    # START
    fs = 14
    # pfun.start_plot(fs=fs, figsize=(20,9))
    pfun.start_plot(fs=fs, figsize=(15,7)) #make smaller for paper
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])

    # PLOT CODE
    vn = 'salt'
    if vn == 'salt':
        pinfo.cmap_dict[vn] = 'Spectral_r'

    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])

    # CREATE THE SECTION
    # create track by hand
    lon = G['lon_rho']
    lat = G['lat_rho']
    zdeep = -205
    #x = np.linspace(1.1, -1, 500) #sill1
    #x_e = np.linspace(0.45, 0.675, 100) #use less points for shorter section?

    #use 2.1 for all to get consistent x-axis scaling
    x_e = np.linspace(0, 2.1, 500) #use less points for shorter section? 1.1 for 5km model, 1.3 for 20km, 2.1 for 80km
    #xcoast = 94.3 #sill2
    y_e = 45 * np.ones(x_e.shape)

    #v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict) #old section function
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = pfun.get_sect(in_dict['fn'], vn, x_e, y_e)

    # COLOR
    # scaled section data
    #sf = pinfo.fac_dict[vn] * v3['sectvarf'] #old
    sf = pinfo.fac_dict[vn] * fld_s #new

    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING
    # # map with section line
    # ax = fig.add_subplot(1, 3, 1)
    # cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
    #         cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    # # fig.colorbar(cs, ax=ax) # It is identical to that of the section
    # #pfun.add_coast(ax)
    # aaf = [0.45, 0.675, 44.95, 45.05] # focus domain
    # #aaf = [-4, 2, 43, 47] # focus domain
    # ax.axis(aaf)
    # pfun.dar(ax)
    # pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    # ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    # ax.set_xlabel('Longitude')
    # ax.set_ylabel('Latitude')
    # # add section track
    # ax.plot(x, y, '-r', linewidth=2)
    # # ax.plot(x[idist0], y[idist0], 'r', markersize=5, markerfacecolor='w',
    # #     markeredgecolor='r', markeredgewidth=2) #old
    # ax.plot(x[0], y[0], 'or', markersize=5, markerfacecolor='w',
    #     markeredgecolor='r', markeredgewidth=2) #new
    # ax.set_xticks([0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675])
    # ax.set_yticks([44.95, 45, 45.05])
    # # ax.set_xticks([-4, -2, 0, 2])
    # # ax.set_yticks([43, 45, 47])

    # section
    ax = fig.add_subplot(1, 1, 1) #no map, only contour section
    # ax.plot(dist, v2['zbot'], 'k', linewidth=2)
    # ax.plot(dist, v2['zeta'], '-k', linewidth=2)
    ax.plot(dist, zbot, '-k', linewidth=2)
    ax.plot(dist, ztop, '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    #ax.invert_xaxis() #not sure if this needs to be changed?
    #ax.set_xticks([0, 2, 4, 6, 8, 10, 12]) #comment out for now
    #ax.vlines(xcoast, -200, 0, linestyles='dashed') #add coast line #don't need since zoomed in
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    # cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:],
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #old
    # cs = ax.pcolormesh(dist_se,zw_se,sf,
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #new
    # cs = ax.contour(dist_se,zw_se,sf,
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #contour
    cs = ax.contourf((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
                        levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], cmap=pinfo.cmap_dict[vn],extend='both') #contour with manual vmax/vmin
    cs = ax.contour((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
                        levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], colors='k') #contour with manual vmax/vmin
    #cs = ax.contourf((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
    #                   vmin=22, vmax=34, levels=25, cmap=pinfo.cmap_dict[vn]) #contour with manual vmax/vmin
    
    ax.clabel(cs, inline=True, fontsize=12)
    #fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance [km]')
    ax.set_ylabel('Z [m]')
    
    # ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn])) #no title for paper plots

    #this is similar to pfun.add info, but does not change the timezone and does not include showing the grid name
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt = T['dt']
    ax.text(.97, .11, dt.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    ax.text(.97, .1, dt.strftime('%H:%M'),
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    #add letters for paper subplots
    gridname=(str(in_dict['fn']).split('/')[-3]).split('_')[0]
    if gridname=='sill5km':
        ax.text(.03, .1, 'A',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    elif gridname=='sill20kmdeep':
        ax.text(.03, .1, 'B',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    elif gridname=='sill80km':
        ax.text(.03, .1, 'C',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_contour_sillslope_eb(in_dict):
    """
    This plots a section (distance, z) with line contours, and makes sure
    that the color limits are identical.

    Change limits so that it is only around the end of the sill, and make equal axis scaling
    
    Uses the new pfun.get_sect() function.
    """
    # START
    fs = 14
    # pfun.start_plot(fs=fs, figsize=(20,9))
    # pfun.start_plot(fs=fs, figsize=(15,7)) #make smaller for paper
    pfun.start_plot(fs=fs, figsize=(20,5)) #make smaller for paper
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])

    # PLOT CODE
    vn = 'salt'
    if vn == 'salt':
        pinfo.cmap_dict[vn] = 'Spectral_r'

    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])

    # CREATE THE SECTION
    # create track by hand
    lon = G['lon_rho']
    lat = G['lat_rho']
    zdeep = -205
    #x = np.linspace(1.1, -1, 500) #sill1
    #x_e = np.linspace(0.45, 0.675, 100) #use less points for shorter section?

    #use 2.1 for all to get consistent x-axis scaling
    x_e = np.linspace(0, 2.1, 500) #use less points for shorter section? 1.1 for 5km model, 1.3 for 20km, 2.1 for 80km
    #xcoast = 94.3 #sill2
    y_e = 45 * np.ones(x_e.shape)

    #v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict) #old section function
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = pfun.get_sect(in_dict['fn'], vn, x_e, y_e)

    # COLOR
    # scaled section data
    #sf = pinfo.fac_dict[vn] * v3['sectvarf'] #old
    sf = pinfo.fac_dict[vn] * fld_s #new

    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING
    # # map with section line
    # ax = fig.add_subplot(1, 3, 1)
    # cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
    #         cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    # # fig.colorbar(cs, ax=ax) # It is identical to that of the section
    # #pfun.add_coast(ax)
    # aaf = [0.45, 0.675, 44.95, 45.05] # focus domain
    # #aaf = [-4, 2, 43, 47] # focus domain
    # ax.axis(aaf)
    # pfun.dar(ax)
    # pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    # ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    # ax.set_xlabel('Longitude')
    # ax.set_ylabel('Latitude')
    # # add section track
    # ax.plot(x, y, '-r', linewidth=2)
    # # ax.plot(x[idist0], y[idist0], 'r', markersize=5, markerfacecolor='w',
    # #     markeredgecolor='r', markeredgewidth=2) #old
    # ax.plot(x[0], y[0], 'or', markersize=5, markerfacecolor='w',
    #     markeredgecolor='r', markeredgewidth=2) #new
    # ax.set_xticks([0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675])
    # ax.set_yticks([44.95, 45, 45.05])
    # # ax.set_xticks([-4, -2, 0, 2])
    # # ax.set_yticks([43, 45, 47])

    # section
    ax = fig.add_subplot(1, 1, 1) #no map, only contour section
    # ax.plot(dist, v2['zbot'], 'k', linewidth=2)
    # ax.plot(dist, v2['zeta'], '-k', linewidth=2)
    ax.plot(dist, zbot, '-k', linewidth=2)
    ax.plot(dist, ztop, '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    #ax.invert_xaxis() #not sure if this needs to be changed?
    #ax.set_xticks([0, 2, 4, 6, 8, 10, 12]) #comment out for now
    #ax.vlines(xcoast, -200, 0, linestyles='dashed') #add coast line #don't need since zoomed in
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    # cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:],
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #old
    # cs = ax.pcolormesh(dist_se,zw_se,sf,
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #new
    # cs = ax.contour(dist_se,zw_se,sf,
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #contour
    cs = ax.contourf((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
                        levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], cmap=pinfo.cmap_dict[vn],extend='both') #contour with manual vmax/vmin
    cs = ax.contour((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
                        levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], colors='k') #contour with manual vmax/vmin
    #cs = ax.contourf((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
    #                   vmin=22, vmax=34, levels=25, cmap=pinfo.cmap_dict[vn]) #contour with manual vmax/vmin
    
    ax.clabel(cs, inline=True, fontsize=12)
    #fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance [km]')
    ax.set_ylabel('Z [m]')
    
    # ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn])) #no title for paper plots

    #this is similar to pfun.add info, but does not change the timezone and does not include showing the grid name
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt = T['dt']
    ax.text(.97, .11, dt.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    ax.text(.97, .1, dt.strftime('%H:%M'),
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    #add letters for paper subplots
    gridname=(str(in_dict['fn']).split('/')[-3]).split('_')[0]
    if gridname=='sill5km':
        ax.text(.03, .1, 'A',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
        ax.set_xlim(44.5,47.5)
    elif gridname=='sill20kmdeep':
        ax.text(.03, .1, 'B',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
        ax.set_xlim(59.5,62.5)
    elif gridname=='sill80km':
        ax.text(.03, .1, 'C',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
        ax.set_xlim(119.5,122.5)

    ax.set_aspect(0.001)
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_contourzoom_eb(in_dict):
    """
    Map and section countours around sill. This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    Uses the new pfun.get_sect() function.
    """
    # START
    # fs = 14
    # pfun.start_plot(fs=fs, figsize=(20,9))
    # fig = plt.figure()
    # fig, axs = plt.subplots(2, 1, figsize=(10,12),gridspec_kw={'height_ratios': [4,1]})
    fig = plt.figure(figsize=(15,12)) #sill2
    gs = fig.add_gridspec(nrows=2,ncols=2, width_ratios=[25,1], height_ratios=[3,1])
    ax0 = fig.add_subplot(gs[0,0])
    ax1 = fig.add_subplot(gs[1,0])
    ax2 = fig.add_subplot(gs[0,1])
    ds = xr.open_dataset(in_dict['fn'])

    # PLOT CODE
    vn = 'salt'
    if vn == 'salt':
        pinfo.cmap_dict[vn] = 'Spectral_r'

    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])

    # CREATE THE SECTION
    # create track by hand
    lon = G['lon_rho']
    lat = G['lat_rho']
    zdeep = -205
    #x = np.linspace(1.1, -1, 500) #sill1
    #x_e = np.linspace(0.45, 0.675, 100) #use less points for shorter section?
    #set east x limit based on length of model 1.1 for 5km model, 1.3 for 20km, 2.1 for 80km
    gridname=(str(in_dict['fn']).split('/')[-3]).split('_')[0]
    if gridname=='sill5km':
        xlonlim=1.1
    elif gridname=='sill10km':
        xlonlim=1.2
    elif gridname=='sill20kmdeep':
        xlonlim=1.3
    elif gridname=='sill40km':
        xlonlim=1.6
    elif gridname=='sill80km':
        xlonlim=2.1

    x_e = np.linspace(0, xlonlim, 500) #use less points for shorter section? 
    aaf = [0, xlonlim, 44.95, 45.05] # focus domain
    #xcoast = 94.3 #sill2
    y_e = 45 * np.ones(x_e.shape)

    #v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict) #old section function
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = pfun.get_sect(in_dict['fn'], vn, x_e, y_e)

    # COLOR
    # scaled section data
    #sf = pinfo.fac_dict[vn] * v3['sectvarf'] #old
    sf = pinfo.fac_dict[vn] * fld_s #new

    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    # if in_dict['auto_vlims']:
    #     pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    pinfo.vlims_dict[vn]=(24,34)
    
    # PLOTTING
    # map with section line
    #ax = fig.add_subplot(1, 3, 1)
    # cs = pfun.add_map_field(ax1, ds, vn, pinfo.vlims_dict,
    #         cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    ax1.contourf(lon,lat,ds[vn][0, -1,:,:].values,
                        levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], cmap=pinfo.cmap_dict[vn],extend='both') #contour with manual vmax/vmin
    ax1.contour(lon,lat,ds[vn][0, -1,:,:].values,
                        levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], colors='k', linewidths=1) #contour with manual vmax/vmin
    # fig.colorbar(cs, ax=ax) # It is identical to that of the section
    #pfun.add_coast(ax)
    #aaf = [-4, 2, 43, 47] # focus domain
    ax1.axis(aaf)
    pfun.dar(ax1)
    ax1.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    # add section track
    ax1.plot(x, y, '-r', linewidth=1)
    # ax.plot(x[idist0], y[idist0], 'r', markersize=5, markerfacecolor='w',
    #     markeredgecolor='r', markeredgewidth=2) #old
    #ax1.plot(x[0], y[0], 'or', markersize=5, markerfacecolor='w',
        #markeredgecolor='r', markeredgewidth=2) #new
    # ax.set_xticks([0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6, 0.625, 0.65, 0.675])
    # ax.set_yticks([44.95, 45, 45.05])
    # ax.set_xticks([-4, -2, 0, 2])
    # ax.set_yticks([43, 45, 47])

    # section
    #ax = fig.add_subplot(1, 3, (2, 3))
    # ax.plot(dist, v2['zbot'], 'k', linewidth=2)
    # ax.plot(dist, v2['zeta'], '-k', linewidth=2)
    ax0.plot(dist, zbot, '-k', linewidth=2)
    ax0.plot(dist, ztop, '-b', linewidth=1)
    ax0.set_xlim(dist.min(), dist.max())
    #ax.invert_xaxis() #not sure if this needs to be changed?
    #ax.set_xticks([0, 2, 4, 6, 8, 10, 12]) #comment out for now
    #ax.vlines(xcoast, -200, 0, linestyles='dashed') #add coast line #don't need since zoomed in
    ax0.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    # cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:],
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #old
    # cs = ax.pcolormesh(dist_se,zw_se,sf,
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #new
    cs = ax0.contourf((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
                        levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], cmap=pinfo.cmap_dict[vn],extend='both') #contour with manual vmax/vmin
    cs2 = ax0.contour((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
                        levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], colors='k') #contour with manual vmax/vmin
    
    fig.colorbar(cs, cax=ax2)
    ax0.set_xlabel('Distance (km)')
    ax0.set_ylabel('Z (m)')
    ax0.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    pfun.add_info(ax0, in_dict['fn'], loc='lower_right')
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_uw_quiver_eb(in_dict):
    """
    This plots a section (distance, z) with line contours, and makes sure
    that the color limits are identical.
    
    Uses the new pfun.get_sect() function.
    """
    # START
    fs = 14
    # pfun.start_plot(fs=fs, figsize=(20,9))
    pfun.start_plot(fs=fs, figsize=(15,7)) #make smaller for paper
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # make u_rho #and add w_rho #maybe can skip this??
    # u_rho = ds.salt * np.nan
    # u = ds.u.values
    # u_rho[:,:,:,1:-1] = (u[:,:,:,1:]+u[:,:,:,:-1])/2
    # ds['u_rho']=u_rho

    # PLOT CODE
    # vn = 'u_rho'
    # if vn == 'u_rho':
    #     pinfo.cmap_dict[vn] = cm.balance
    #     pinfo.fac_dict[vn] = 1
    vn = 'u'
    if vn == 'u':
        pinfo.cmap_dict[vn] = cm.balance

    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])

    # CREATE THE SECTION
    # create track by hand
    lon = G['lon_rho']
    lat = G['lat_rho']
    zdeep = -205
    #x = np.linspace(1.1, -1, 500) #sill1
    #x_e = np.linspace(0.45, 0.675, 100) #use less points for shorter section?

    #use 2.1 for all to get consistent x-axis scaling
    x_e = np.linspace(0, 2.1, 500) #use less points for shorter section? 1.1 for 5km model, 1.3 for 20km, 2.1 for 80km
    #xcoast = 94.3 #sill2
    y_e = 45 * np.ones(x_e.shape)

    #v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict) #old section function
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = pfun.get_sect(in_dict['fn'], vn, x_e, y_e)

    #get another section with w
    vn2 = 'w'
    x2, y2, dist2, dist_e2, zbot2, ztop2, dist_se2, zw_se2, fld_s2, lon2, lat2 = pfun.get_sect(in_dict['fn'], vn2, x_e, y_e)

    # COLOR
    # scaled section data
    #sf = pinfo.fac_dict[vn] * v3['sectvarf'] #old
    sf_u = pinfo.fac_dict[vn] * fld_s #new
    sf_w = pinfo.fac_dict[vn2] * fld_s2

    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf_u)
    
    # PLOTTING
    # section
    ax = fig.add_subplot(1, 1, 1) #no map, only contour section
    # ax.plot(dist, v2['zbot'], 'k', linewidth=2)
    # ax.plot(dist, v2['zeta'], '-k', linewidth=2)
    ax.plot(dist, zbot, '-k', linewidth=2)
    ax.plot(dist, ztop, '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    #ax.invert_xaxis() #not sure if this needs to be changed?
    #ax.set_xticks([0, 2, 4, 6, 8, 10, 12]) #comment out for now
    #ax.vlines(xcoast, -200, 0, linestyles='dashed') #add coast line #don't need since zoomed in
    ax.set_ylim(zdeep, 5)
    # plot section
    # svlims = pinfo.vlims_dict[vn]
    # cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:],
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #old
    # cs = ax.pcolormesh(dist_se,zw_se,sf,
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #new
    # cs = ax.contour(dist_se,zw_se,sf,
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #contour
    # cs = ax.contourf((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
    #                     levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], cmap=pinfo.cmap_dict[vn],extend='both') #contour with manual vmax/vmin
    # cs = ax.contour((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
    #                     levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], colors='k') #contour with manual vmax/vmin
    #cs = ax.contourf((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
    #                   vmin=22, vmax=34, levels=25, cmap=pinfo.cmap_dict[vn]) #contour with manual vmax/vmin
    
    #pcolormesh u plot
    # vmin = -0.1
    # vmax = 0.1
    # cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:], vmin=vmin, vmax=vmax, cmap=pinfo.cmap_dict[vn])
    # cs = ax.pcolormesh(dist_se, zw_se, sf, vmin=vmin, vmax=vmax, cmap=pinfo.cmap_dict[vn])
    #contour u plot
    dist_plot = 1000*(dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4 #the plotting distance in m instead of km (need to do this for using quiver more easily)
    z_plot = (zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4
    cs = ax.contourf(dist_plot,z_plot,sf_u,
                        levels=[-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2], cmap=pinfo.cmap_dict[vn],extend='both') #contour with manual vmax/vmin

    # ax.clabel(cs, inline=True, fontsize=12)#, color='tab:gray')
    fig.colorbar(cs,ax=ax,label=r'$u$ [m/s]',location='bottom') #turn this on to make a plot just for the colorbar

    #add quiver plot
    Q = ax.quiver(dist_plot[::2,15::20], z_plot[::2,15::20], sf_u[::2,15::20], sf_w[::2,15::20], #NEED TO SUBSAMPLE THE ARROWS!!
              angles='uv', color='k') # scale=2, scale_units='inches', units='inches', #NEED TO SET THE SCALE SO IT WILL BE THE SAME FOR ALL PLOTS
    qk = ax.quiverkey(Q, 0.9, 0.1, 0.2, r'$0.2$ m/s', labelpos='E', coordinates='figure')

    #fig.colorbar(cs, ax=ax)
    #fix the units on the x axis since we plotted in m
    ax.set_xticks([0,20e3,40e3,60e3,80e3,100e3,120e3,140e3,160e3],['0','20','40','60','80','100','120','140','160'])
    ax.set_xlabel('Distance [km]')
    ax.set_ylabel('Z [m]')
    
    # ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn])) #no title for paper plots

    #this is similar to pfun.add info, but does not change the timezone and does not include showing the grid name
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt = T['dt']
    ax.text(.97, .11, dt.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    ax.text(.97, .1, dt.strftime('%H:%M'),
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    #add letters for paper subplots
    gridname=(str(in_dict['fn']).split('/')[-3]).split('_')[0]
    if gridname=='sill5km':
        ax.text(.03, .1, 'A',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    elif gridname=='sill20kmdeep':
        ax.text(.03, .1, 'B',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    elif gridname=='sill80km':
        ax.text(.03, .1, 'C',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_contour_u_eb(in_dict):
    """
    This plots a section (distance, z) with line contours, and makes sure
    that the color limits are identical.
    
    Uses the new pfun.get_sect() function.
    """
    # START
    fs = 14
    # pfun.start_plot(fs=fs, figsize=(20,9))
    pfun.start_plot(fs=fs, figsize=(15,7)) #make smaller for paper
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # make u_rho #and add w_rho #maybe can skip this??
    # u_rho = ds.salt * np.nan
    # u = ds.u.values
    # u_rho[:,:,:,1:-1] = (u[:,:,:,1:]+u[:,:,:,:-1])/2
    # ds['u_rho']=u_rho

    # PLOT CODE
    # vn = 'u_rho'
    # if vn == 'u_rho':
    #     pinfo.cmap_dict[vn] = cm.balance
    #     pinfo.fac_dict[vn] = 1
    vn = 'u'
    if vn == 'u':
        # pinfo.cmap_dict[vn] = cm.balance
        pinfo.cmap_dict[vn] = 'RdBu_r'

    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])

    # CREATE THE SECTION
    # create track by hand
    lon = G['lon_rho']
    lat = G['lat_rho']
    zdeep = -205
    #x = np.linspace(1.1, -1, 500) #sill1
    #x_e = np.linspace(0.45, 0.675, 100) #use less points for shorter section?

    #use 2.1 for all to get consistent x-axis scaling
    x_e = np.linspace(0, 2.1, 500) #use less points for shorter section? 1.1 for 5km model, 1.3 for 20km, 2.1 for 80km
    #xcoast = 94.3 #sill2
    y_e = 45 * np.ones(x_e.shape)

    #v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict) #old section function
    x, y, dist, dist_e, zbot, ztop, dist_se, zw_se, fld_s, lon, lat = pfun.get_sect(in_dict['fn'], vn, x_e, y_e)

    # COLOR
    # scaled section data
    #sf = pinfo.fac_dict[vn] * v3['sectvarf'] #old
    sf_u = pinfo.fac_dict[vn] * fld_s #new

    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf_u)
    
    # PLOTTING
    # section
    ax = fig.add_subplot(1, 1, 1) #no map, only contour section
    # ax.plot(dist, v2['zbot'], 'k', linewidth=2)
    # ax.plot(dist, v2['zeta'], '-k', linewidth=2)
    ax.plot(dist, zbot, '-k', linewidth=2)
    ax.plot(dist, ztop, '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    #ax.invert_xaxis() #not sure if this needs to be changed?
    #ax.set_xticks([0, 2, 4, 6, 8, 10, 12]) #comment out for now
    #ax.vlines(xcoast, -200, 0, linestyles='dashed') #add coast line #don't need since zoomed in
    ax.set_ylim(zdeep, 5)
    # plot section
    # svlims = pinfo.vlims_dict[vn]
    # cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:],
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #old
    # cs = ax.pcolormesh(dist_se,zw_se,sf,
    #                    vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn]) #new
    # cs = ax.contourf((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
    #                     levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], cmap=pinfo.cmap_dict[vn],extend='both') #contour with manual vmax/vmin
    # cs = ax.contour((dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4,(zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4,sf,
    #                     levels=[24,24.5,25,25.5,26,26.5,27,27.5,28,28.5,29,29.5,30,30.5,31,31.5,32,32.5,33,33.5,34], colors='k') #contour with manual vmax/vmin
    
    #pcolormesh u plot
    # vmin = -0.1
    # vmax = 0.1
    # cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:], vmin=vmin, vmax=vmax, cmap=pinfo.cmap_dict[vn])
    # cs = ax.pcolormesh(dist_se, zw_se, sf, vmin=vmin, vmax=vmax, cmap=pinfo.cmap_dict[vn])
    #contour u plot
    dist_plot = (dist_se[:-1,:-1]+dist_se[1:,:-1]+dist_se[:-1,1:]+dist_se[1:,1:])/4
    z_plot = (zw_se[:-1,:-1]+zw_se[1:,:-1]+zw_se[:-1,1:]+zw_se[1:,1:])/4
    csf = ax.contourf(dist_plot,z_plot,sf_u,
                        levels=[-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2], cmap=pinfo.cmap_dict[vn],extend='both') #contour with manual vmax/vmin
    cs = ax.contour(dist_plot,z_plot,sf_u,
                        levels=[-0.2,-0.18,-0.16,-0.14,-0.12,-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2], colors='k',linewidths=0.3, linestyles='solid',negative_linestyles='solid') #contour with manual vmax/vmin
    # ax.clabel(cs, inline=True,levels=[-0.2,-0.16,-0.12,-0.08,-0.04,0,0.04,0.08,0.12,0.16,0.2], fontsize=8)
    fig.colorbar(csf,ax=ax,label=r'$u$ [m/s]',location='bottom',fraction=0.05,shrink=0.6) #turn this on to make a plot just for the colorbar
    ax.set_xlabel('Distance [km]')
    ax.set_ylabel('Z [m]')
    
    # ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn])) #no title for paper plots

    #this is similar to pfun.add info, but does not change the timezone and does not include showing the grid name
    T = zrfun.get_basic_info(in_dict['fn'], only_T=True)
    dt = T['dt']
    ax.text(.97, .11, dt.strftime('%Y-%m-%d'),
            horizontalalignment='right' , verticalalignment='bottom',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    ax.text(.97, .1, dt.strftime('%H:%M'),
            horizontalalignment='right', verticalalignment='top',
            transform=ax.transAxes, fontsize=fs,
            bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    #add letters for paper subplots
    gridname=(str(in_dict['fn']).split('/')[-3]).split('_')[0]
    if gridname=='sill5km':
        ax.text(.03, .1, 'A',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    elif gridname=='sill20kmdeep':
        ax.text(.03, .1, 'B',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    elif gridname=='sill80km':
        ax.text(.03, .1, 'C',
            horizontalalignment='left', verticalalignment='top',
            transform=ax.transAxes, fontsize=24, fontweight='bold')#,
            #bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_u_eb(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'u_rho'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    u_rho = ds.salt * np.nan
    u = ds.u.values
    u_rho[:,:,:,1:-1] = (u[:,:,:,1:]+u[:,:,:,:-1])/2
    ds['u_rho']=u_rho
    # CREATE THE SECTION
    # create track by hand
    lon = G['lon_rho']
    lat = G['lat_rho']
    zdeep = -205
    #x = np.linspace(1.1, -1, 500) #sill1
    x = np.linspace(1.2, -1, 500) #sill2
    xcoast = 94.3 #sill2
    y = 45 * np.ones(x.shape)
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    fac_dict =  {'u_rho': 1}
    sf = fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    cmap = copy.copy(cm.balance)
    cmap.set_bad('lightgray')
    vmin = -0.5
    vmax = 0.5
    vlims_dict = {'u_rho': (vmin, vmax)}
    
    # PLOTTING
    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    cs = pfun.add_map_field(ax, ds, vn, vlims_dict,
            cmap=cmap, fac=fac_dict[vn], do_mask_edges=True)
    # fig.colorbar(cs, ax=ax) # It is identical to that of the section
    #pfun.add_coast(ax)
    aaf = [-4, 4, 43, 47] # focus domain
    ax.axis(aaf)
    pfun.dar(ax)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface u velocity [m/s]')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'r', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    ax.set_xticks([-4, 0, 4])
    ax.set_yticks([43, 45, 47])
    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], 'k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-k', linewidth=2)
    ax.set_xlim(dist.min(), dist.max())
    ax.invert_xaxis()
    ax.set_xticks([0, 40, 80, 120, 160, 200])
    ax.set_ylim(zdeep, 5)
    ax.vlines(xcoast, -200, 0, linestyles='dashed', colors='k') #add coast line
    # plot section
    #svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:],
                       vmin=vmin, vmax=vmax, cmap=cmap)
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section u velocity [m/s]')
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_sect_v_eb(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'v_rho'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    v_rho = ds.salt * np.nan
    v = ds.v.values
    v_rho[:,:,1:-1,:] = (v[:,:,1:,:]+v[:,:,:-1,:])/2
    ds['v_rho']=v_rho
    # CREATE THE SECTION
    # create track by hand
    lon = G['lon_rho']
    lat = G['lat_rho']
    zdeep = -205
    #x = np.linspace(1.1, -1, 500) #sill1
    x = np.linspace(1.2, -1, 500) #sill2
    xcoast = 94.3 #sill2
    y = 45 * np.ones(x.shape)
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    fac_dict =  {'v_rho': 1}
    sf = fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    cmap = copy.copy(cm.balance)
    cmap.set_bad('lightgray')
    vmin = -0.5
    vmax = 0.5
    vlims_dict = {'v_rho': (vmin, vmax)}
    
    # PLOTTING
    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    cs = pfun.add_map_field(ax, ds, vn, vlims_dict,
            cmap=cmap, fac=fac_dict[vn], do_mask_edges=True)
    # fig.colorbar(cs, ax=ax) # It is identical to that of the section
    #pfun.add_coast(ax)
    aaf = [-4, 4, 43, 47] # focus domain
    ax.axis(aaf)
    pfun.dar(ax)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface v velocity [m/s]')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'r', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    ax.set_xticks([-4, 0, 4])
    ax.set_yticks([43, 45, 47])
    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], 'k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-k', linewidth=2)
    ax.set_xlim(dist.min(), dist.max())
    ax.invert_xaxis()
    ax.set_xticks([0, 40, 80, 120, 160, 200])
    ax.set_ylim(zdeep, 5)
    ax.vlines(xcoast, -200, 0, linestyles='dashed', colors='k') #add coast line
    # plot section
    #svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'][1:-1,:], v3['zrf'][1:-1,:], sf[1:-1,:],
                       vmin=vmin, vmax=vmax, cmap=cmap)
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section v velocity [m/s]')
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_vort_eb(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    aa1 = [-2, 0, 44.5, 46.5]
    #aa2 = [-4, 4, 43, 47] #sill1
    aa2 = [-4, 2, 43, 47] #sill2
    aa3 = [-0.2, 1.2, 44.9, 45.1]
    # find aspect ratio of the map
    # AR is the aspect ratio of the map: Vertical/Horizontal
    # AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    # fs = 14
    # hgt = 10
    #pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))

    
    # create fields
    u = ds.u[0,-1,:,:].values
    v = ds.v[0,-1,:,:].values
    dx = 1/ds.pm.values
    dy = 1/ds.pn.values
    # vort is on the psi grid (plot with lon_rho, lat_rho)
    vort = np.diff(v,axis=1)/dx[1:,1:] - np.diff(u,axis=0)/dy[1:,1:]
    # dive is on the trimmed rho grid
    # dive = np.diff(u[1:-1,:], axis=1)/dx[1:-1,1:-1] + np.diff(v[:,1:-1],axis=0)/dy[1:-1,1:-1]

    
    # set color limits
    vv = 2*np.nanstd(vort)
    if vv<=0.0002:
        vv=0.0002

    #fig = plt.figure(figsize=(14,8)) #sill1
    # gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[17,10,1], height_ratios=[3,1]) #sill1
    fig = plt.figure(figsize=(12,8)) #sill2
    gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[13,10,1], height_ratios=[3,1]) #sill2
    cmap = copy.copy(plt.cm.PiYG_r)
    cmap.set_bad('lightgray')


    # PLOT CODE
    if in_dict['auto_vlims']: #MUST USE -avl True SINCE vort IS NOT INCLUDED IN pinfo
        pinfo.vlims_dict['vort'] = (-vv, vv)
        pinfo.vlims_dict['dive'] = (-vv, vv)

    vmin = pinfo.vlims_dict['vort'][0]
    vmax = pinfo.vlims_dict['vort'][1]
    
    ax1 = fig.add_subplot(gs[0,1]) 
    cs1 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, vort, cmap=cmap, vmin = vmin, vmax = vmax)
    ax1.set_title('Plume focus', fontsize=12)
    #fig.colorbar(cs1)
    ax1.axis(aa1)
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')

    ax2 = fig.add_subplot(gs[0,0])
    cs2 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, vort, cmap=cmap, vmin = vmin, vmax = vmax)
    ax2.set_title('Full model', fontsize=12)
    #fig.colorbar(cs2)
    ax2.axis(aa2)
    pfun.dar(ax2)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    pfun.add_info(ax2, in_dict['fn'])
    # pfun.add_coast(ax2)
    # pfun.add_bathy_contours(ax, ds, txt=True)

    ax3 = fig.add_subplot(gs[1,0:2])
    cs3 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, vort, cmap=cmap, vmin = vmin, vmax = vmax)
    ax3.set_title('Estuary focus', fontsize=12)
    ax3.axis(aa3)
    pfun.dar(ax3)
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')

    ax4 = fig.add_subplot(gs[:,2])
    fig.colorbar(cs3, cax=ax4)
    plt.suptitle('Surface Vorticity $[s^{-1}]$', fontsize=16)
    #plt.tight_layout()

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_dive_eb(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    aa1 = [-2, 0, 44.5, 46.5]
    # aa2 = [-4, 4, 43, 47] #sill1
    aa2 = [-4, 2, 43, 47] #sill2
    aa3 = [-0.2, 1.2, 44.9, 45.1]
    # find aspect ratio of the map
    # AR is the aspect ratio of the map: Vertical/Horizontal
    # AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    # fs = 14
    # hgt = 10
    #pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))

    
    # create fields
    u = ds.u[0,-1,:,:].values
    v = ds.v[0,-1,:,:].values
    dx = 1/ds.pm.values
    dy = 1/ds.pn.values
    # vort is on the psi grid (plot with lon_rho, lat_rho)
    # vort = np.diff(v,axis=1)/dx[1:,1:] - np.diff(u,axis=0)/dy[1:,1:]
    # dive is on the trimmed rho grid
    dive = np.diff(u[1:-1,:], axis=1)/dx[1:-1,1:-1] + np.diff(v[:,1:-1],axis=0)/dy[1:-1,1:-1]

    
    # set color limits
    vv = 2*np.nanstd(dive)
    if vv<=0.0002:
        vv=0.0002

    #fig = plt.figure(figsize=(14,8)) #sill1
    #gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[17,10,1], height_ratios=[3,1]) #sill1
    fig = plt.figure(figsize=(12,8)) #sill2
    gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[13,10,1], height_ratios=[3,1]) #sill2
    cmap = copy.copy(plt.cm.BrBG_r)
    cmap.set_bad('lightgray')


    # PLOT CODE
    if in_dict['auto_vlims']: #MUST USE -avl True SINCE vort IS NOT INCLUDED IN pinfo
        pinfo.vlims_dict['dive'] = (-vv, vv)

    vmin = pinfo.vlims_dict['dive'][0]
    vmax = pinfo.vlims_dict['dive'][1]

    plon, plat = pfun.get_plon_plat(ds.lon_rho[1:-1,1:-1].values, ds.lat_rho[1:-1,1:-1].values)
    
    ax1 = fig.add_subplot(gs[0,1]) 
    cs1 = plt.pcolormesh(plon, plat, dive, cmap=cmap, vmin = vmin, vmax = vmax)
    ax1.set_title('Plume focus', fontsize=12)
    #fig.colorbar(cs1)
    ax1.axis(aa1)
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')

    ax2 = fig.add_subplot(gs[0,0])
    cs2 = plt.pcolormesh(plon, plat, dive, cmap=cmap, vmin = vmin, vmax = vmax)
    ax2.set_title('Full model', fontsize=12)
    #fig.colorbar(cs2)
    ax2.axis(aa2)
    pfun.dar(ax2)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    pfun.add_info(ax2, in_dict['fn'])
    # pfun.add_coast(ax2)
    # pfun.add_bathy_contours(ax, ds, txt=True)

    ax3 = fig.add_subplot(gs[1,0:2])
    cs3 = plt.pcolormesh(plon, plat, dive, cmap=cmap, vmin = vmin, vmax = vmax)
    ax3.set_title('Estuary focus', fontsize=12)
    ax3.axis(aa3)
    pfun.dar(ax3)
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')

    ax4 = fig.add_subplot(gs[:,2])
    fig.colorbar(cs3, cax=ax4)
    plt.suptitle('Surface Divergence $[s^{-1}]$', fontsize=16)
    #plt.tight_layout()

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_saltmap_eb(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    aa1 = [-2, 0, 44.5, 46.5]
    # aa2 = [-4, 4, 43, 47] #sill1
    aa2 = [-4, 2, 43, 47] #sill2
    aa3 = [-0.2, 1.3, 44.9, 45.1]
    # find aspect ratio of the map
    # AR is the aspect ratio of the map: Vertical/Horizontal
    # AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    # fs = 14
    # hgt = 10
    #pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))

    
    # create fields
    salt = ds.salt[0,-1,:,:].values
    # u = ds.u[0,-1,:,:].values
    # v = ds.v[0,-1,:,:].values
    # dx = 1/ds.pm.values
    # dy = 1/ds.pn.values
    # # vort is on the psi grid (plot with lon_rho, lat_rho)
    # vort = np.diff(v,axis=1)/dx[1:,1:] - np.diff(u,axis=0)/dy[1:,1:]
    # dive is on the trimmed rho grid
    # dive = np.diff(u[1:-1,:], axis=1)/dx[1:-1,1:-1] + np.diff(v[:,1:-1],axis=0)/dy[1:-1,1:-1]

    # # set color limits
    # vv = 2*np.nanstd(vort)
    # if vv<=0.00001:
    #     vv=0.0002

    #fig = plt.figure(figsize=(14,8)) #sill1
    #gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[17,10,1], height_ratios=[3,1]) #sill1
    fig = plt.figure(figsize=(12,8)) #sill2
    gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[13,10,1], height_ratios=[3,1]) #sill2
    cmap = copy.copy(cm.haline)
    cmap.set_bad('lightgray')


    # # PLOT CODE
    # if in_dict['auto_vlims']: #MUST USE -avl True SINCE vort IS NOT INCLUDED IN pinfo
    #     pinfo.vlims_dict['vort'] = (-vv, vv)
    #     pinfo.vlims_dict['dive'] = (-vv, vv)

    # vmin = pinfo.vlims_dict['vort'][0]
    # vmax = pinfo.vlims_dict['vort'][1]

    vmin=20
    vmax=30
    
    ax1 = fig.add_subplot(gs[0,1]) 
    cs1 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, salt, cmap=cmap, vmin = vmin, vmax = vmax)
    ax1.set_title('Plume focus', fontsize=12)
    #fig.colorbar(cs1)
    ax1.axis(aa1)
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')

    ax2 = fig.add_subplot(gs[0,0])
    cs2 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, salt, cmap=cmap, vmin = vmin, vmax = vmax)
    ax2.set_title('Full model', fontsize=12)
    #fig.colorbar(cs2)
    ax2.axis(aa2)
    pfun.dar(ax2)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    pfun.add_info(ax2, in_dict['fn'])
    # pfun.add_coast(ax2)
    # pfun.add_bathy_contours(ax, ds, txt=True)

    ax3 = fig.add_subplot(gs[1,0:2])
    cs3 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, salt, cmap=cmap, vmin = vmin, vmax = vmax)
    ax3.set_title('Estuary focus', fontsize=12)
    ax3.axis(aa3)
    pfun.dar(ax3)
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')

    ax4 = fig.add_subplot(gs[:,2])
    fig.colorbar(cs3, cax=ax4)
    plt.suptitle('Surface Salinity [psu]', fontsize=16)
    #plt.tight_layout()

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_umap_eb(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    aa1 = [-2, 0, 44.5, 46.5]
    # aa2 = [-4, 4, 43, 47] #sill1
    aa2 = [-4, 2, 43, 47] #sill2
    aa3 = [-0.2, 1.2, 44.9, 45.1]
    
    # create fields
    u = ds.u[0,-1,:,:].values

    # set color limits
    # vv = 2*np.nanstd(u)
    # if vv<=0.001:
    #     vv=0.25
    vv=0.25

    #fig = plt.figure(figsize=(14,8)) #sill1
    #gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[17,10,1], height_ratios=[3,1]) #sill1
    fig = plt.figure(figsize=(12,8)) #sill2
    gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[13,10,1], height_ratios=[3,1]) #sill2
    cmap = copy.copy(cm.balance)
    cmap.set_bad('lightgray')


    # # PLOT CODE
    if in_dict['auto_vlims']: #MUST USE -avl True SINCE u IS NOT INCLUDED IN pinfo
        pinfo.vlims_dict['u'] = (-vv, vv)

    vmin = pinfo.vlims_dict['u'][0]
    vmax = pinfo.vlims_dict['u'][1]
    
    ax1 = fig.add_subplot(gs[0,1]) 
    cs1 = plt.pcolormesh(ds.lon_u.values, ds.lat_u.values, u, cmap=cmap, vmin = vmin, vmax = vmax)
    ax1.set_title('Plume focus', fontsize=12)
    #fig.colorbar(cs1)
    ax1.axis(aa1)
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')

    ax2 = fig.add_subplot(gs[0,0])
    cs2 = plt.pcolormesh(ds.lon_u.values, ds.lat_u.values, u, cmap=cmap, vmin = vmin, vmax = vmax)
    ax2.set_title('Full model', fontsize=12)
    #fig.colorbar(cs2)
    ax2.axis(aa2)
    pfun.dar(ax2)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    pfun.add_info(ax2, in_dict['fn'])
    # pfun.add_coast(ax2)
    # pfun.add_bathy_contours(ax, ds, txt=True)

    ax3 = fig.add_subplot(gs[1,0:2])
    cs3 = plt.pcolormesh(ds.lon_u.values, ds.lat_u.values, u, cmap=cmap, vmin = vmin, vmax = vmax)
    ax3.set_title('Estuary focus', fontsize=12)
    ax3.axis(aa3)
    pfun.dar(ax3)
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')

    ax4 = fig.add_subplot(gs[:,2])
    fig.colorbar(cs3, cax=ax4)
    plt.suptitle('Surface u velocity [m/s]', fontsize=16)
    #plt.tight_layout()

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_vmap_eb(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    aa1 = [-2, 0, 44.5, 46.5]
    # aa2 = [-4, 4, 43, 47] #sill1
    aa2 = [-4, 2, 43, 47] #sill2
    aa3 = [-0.2, 1.2, 44.9, 45.1]
    
    # create fields
    v = ds.v[0,-1,:,:].values

    # set color limits
    vv=0.25

    #fig = plt.figure(figsize=(14,8)) #sill1
    #gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[17,10,1], height_ratios=[3,1]) #sill1
    fig = plt.figure(figsize=(12,8)) #sill2
    gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[13,10,1], height_ratios=[3,1]) #sill2
    cmap = copy.copy(cm.balance)
    cmap.set_bad('lightgray')


    # # PLOT CODE
    if in_dict['auto_vlims']: #MUST USE -avl True SINCE u IS NOT INCLUDED IN pinfo
        pinfo.vlims_dict['v'] = (-vv, vv)

    vmin = pinfo.vlims_dict['v'][0]
    vmax = pinfo.vlims_dict['v'][1]
    
    ax1 = fig.add_subplot(gs[0,1]) 
    cs1 = plt.pcolormesh(ds.lon_v.values, ds.lat_v.values, v, cmap=cmap, vmin = vmin, vmax = vmax)
    ax1.set_title('Plume focus', fontsize=12)
    #fig.colorbar(cs1)
    ax1.axis(aa1)
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')

    ax2 = fig.add_subplot(gs[0,0])
    cs2 = plt.pcolormesh(ds.lon_v.values, ds.lat_v.values, v, cmap=cmap, vmin = vmin, vmax = vmax)
    ax2.set_title('Full model', fontsize=12)
    #fig.colorbar(cs2)
    ax2.axis(aa2)
    pfun.dar(ax2)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    pfun.add_info(ax2, in_dict['fn'])
    # pfun.add_coast(ax2)
    # pfun.add_bathy_contours(ax, ds, txt=True)

    ax3 = fig.add_subplot(gs[1,0:2])
    cs3 = plt.pcolormesh(ds.lon_v.values, ds.lat_v.values, v, cmap=cmap, vmin = vmin, vmax = vmax)
    ax3.set_title('Estuary focus', fontsize=12)
    ax3.axis(aa3)
    pfun.dar(ax3)
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')

    ax4 = fig.add_subplot(gs[:,2])
    fig.colorbar(cs3, cax=ax4)
    plt.suptitle('Surface v velocity [m/s]', fontsize=16)
    #plt.tight_layout()

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_spdmap_eb(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    aa1 = [-2, 0, 44.5, 46.5]
    # aa2 = [-4, 4, 43, 47] #sill1
    aa2 = [-4, 2, 43, 47] #sill2
    aa3 = [-0.2, 1.2, 44.9, 45.1]
    
    # create fields
    u_rho = ds.zeta[0,:,:].values * np.nan
    v_rho = ds.zeta[0,:,:].values * np.nan
    u_rho[:,1:-1] = (ds.u[0,-1,:,1:].values + ds.u[0,-1,:,:-1].values)/2
    v_rho[1:-1,:] = (ds.v[0,-1,1:,:].values + ds.v[0,-1,:-1,:].values)/2
    spd = np.sqrt(u_rho**2+v_rho**2)

    # set color limits
    vv=0.5

    #fig = plt.figure(figsize=(14,8)) #sill1
    #gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[17,10,1], height_ratios=[3,1]) #sill1
    fig = plt.figure(figsize=(12,8)) #sill2
    gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[13,10,1], height_ratios=[3,1]) #sill2
    cmap = copy.copy(plt.cm.YlOrBr)
    cmap.set_bad('lightgray')


    # # PLOT CODE
    if in_dict['auto_vlims']: #MUST USE -avl True SINCE u IS NOT INCLUDED IN pinfo
        pinfo.vlims_dict['spd'] = (0, vv)

    vmin = pinfo.vlims_dict['spd'][0]
    vmax = pinfo.vlims_dict['spd'][1]
    
    ax1 = fig.add_subplot(gs[0,1]) 
    cs1 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, spd, cmap=cmap, vmin = vmin, vmax = vmax)
    ax1.quiver(ds.lon_rho[::10,::10].values, ds.lat_rho[::10,::10].values, u_rho[::10,::10], v_rho[::10,::10], color='b', scale=2, scale_units='inches', units='inches', width=0.012)
    ax1.set_title('Plume focus', fontsize=12)
    #fig.colorbar(cs1)
    ax1.axis(aa1)
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')

    ax2 = fig.add_subplot(gs[0,0])
    cs2 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, spd, cmap=cmap, vmin = vmin, vmax = vmax)
    q=ax2.quiver(ds.lon_rho[::15,::15].values, ds.lat_rho[::15,::15].values, u_rho[::15,::15], v_rho[::15,::15], color='b', scale=2, scale_units='inches', units='inches', width=0.012)
    plt.quiverkey(q,0.95,0.9,0.5,'0.5 m/s', angle=45)
    ax2.set_title('Full model', fontsize=12)
    #fig.colorbar(cs2)
    ax2.axis(aa2)
    pfun.dar(ax2)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    pfun.add_info(ax2, in_dict['fn'])
    # pfun.add_coast(ax2)
    # pfun.add_bathy_contours(ax, ds, txt=True)

    ax3 = fig.add_subplot(gs[1,0:2])
    cs3 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, spd, cmap=cmap, vmin = vmin, vmax = vmax)
    ax3.quiver(ds.lon_rho[::4,::4].values, ds.lat_rho[::4,::4].values, u_rho[::4,::4], v_rho[::4,::4], color='b', scale=2, scale_units='inches', units='inches', width=0.012)
    ax3.set_title('Estuary focus', fontsize=12)
    ax3.axis(aa3)
    pfun.dar(ax3)
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')

    ax4 = fig.add_subplot(gs[:,2])
    fig.colorbar(cs3, cax=ax4)
    plt.suptitle('Surface speed [m/s]', fontsize=16)
    #plt.tight_layout()

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_spdbarmap_eb(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    aa1 = [-2, 0, 44.5, 46.5]
    # aa2 = [-4, 4, 43, 47] #sill1
    aa2 = [-4, 2, 43, 47] #sill2
    aa3 = [-0.2, 1.2, 44.9, 45.1]
    
    # create fields
    ubar_rho = ds.zeta[0,:,:].values * np.nan
    vbar_rho = ds.zeta[0,:,:].values * np.nan
    ubar_rho[:,1:-1] = (ds.ubar[0,:,1:].values + ds.ubar[0,:,:-1].values)/2
    vbar_rho[1:-1,:] = (ds.vbar[0,1:,:].values + ds.vbar[0,:-1,:].values)/2
    spdbar = np.sqrt(ubar_rho**2+vbar_rho**2)

    # set color limits
    vv=0.25

    #fig = plt.figure(figsize=(14,8)) #sill1
    #gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[17,10,1], height_ratios=[3,1]) #sill1
    fig = plt.figure(figsize=(12,8)) #sill2
    gs = fig.add_gridspec(nrows=2,ncols=3, width_ratios=[13,10,1], height_ratios=[3,1]) #sill2
    cmap = copy.copy(plt.cm.YlOrBr)
    cmap.set_bad('lightgray')


    # # PLOT CODE
    if in_dict['auto_vlims']: #MUST USE -avl True SINCE u IS NOT INCLUDED IN pinfo
        pinfo.vlims_dict['spdbar'] = (0, vv)

    vmin = pinfo.vlims_dict['spdbar'][0]
    vmax = pinfo.vlims_dict['spdbar'][1]
    
    ax1 = fig.add_subplot(gs[0,1]) 
    cs1 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, spdbar, cmap=cmap, vmin = vmin, vmax = vmax)
    ax1.quiver(ds.lon_rho[::10,::10].values, ds.lat_rho[::10,::10].values, ubar_rho[::10,::10], vbar_rho[::10,::10], color='b', scale=1, scale_units='inches', units='inches', width=0.012)
    ax1.set_title('Plume focus', fontsize=12)
    #fig.colorbar(cs1)
    ax1.axis(aa1)
    pfun.dar(ax1)
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')

    ax2 = fig.add_subplot(gs[0,0])
    cs2 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, spdbar, cmap=cmap, vmin = vmin, vmax = vmax)
    q=ax2.quiver(ds.lon_rho[::15,::15].values, ds.lat_rho[::15,::15].values, ubar_rho[::15,::15], vbar_rho[::15,::15], color='b', scale=1, scale_units='inches', units='inches', width=0.012)
    plt.quiverkey(q,0.95,0.95,0.25,'0.25 m/s', angle=45)
    ax2.set_title('Full model', fontsize=12)
    #fig.colorbar(cs2)
    ax2.axis(aa2)
    pfun.dar(ax2)
    ax2.set_xlabel('Longitude')
    ax2.set_ylabel('Latitude')
    pfun.add_info(ax2, in_dict['fn'])
    # pfun.add_coast(ax2)
    # pfun.add_bathy_contours(ax, ds, txt=True)

    ax3 = fig.add_subplot(gs[1,0:2])
    cs3 = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, spdbar, cmap=cmap, vmin = vmin, vmax = vmax)
    ax3.quiver(ds.lon_rho[::4,::4].values, ds.lat_rho[::4,::4].values, ubar_rho[::4,::4], vbar_rho[::4,::4], color='b', scale=1, scale_units='inches', units='inches', width=0.012)
    ax3.set_title('Estuary focus', fontsize=12)
    ax3.axis(aa3)
    pfun.dar(ax3)
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')

    ax4 = fig.add_subplot(gs[:,2])
    fig.colorbar(cs3, cax=ax4, label=r'$\sqrt{\bar{u}^2+\bar{v}^2}$')
    plt.suptitle('Speed [m/s]', fontsize=16)
    #plt.tight_layout()

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_debug_eb(in_dict):
    # Focused on debugging
    vn_list = ['u', 'v', 'zeta']
    do_wetdry = False
    
    # START
    fs = 10
    pfun.start_plot(fs=fs, figsize=(8*len(vn_list),10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    ii = 1
    for vn in vn_list:
        if 'lon_rho' in ds[vn].coords:
            tag = 'rho'
        if 'lon_u' in ds[vn].coords:
            tag = 'u'
        if 'lon_v' in ds[vn].coords:
            tag = 'v'
        x = ds['lon_'+tag].values
        y = ds['lat_'+tag].values
        px, py = pfun.get_plon_plat(x,y)
        if vn in ['u', 'v']:
            v = ds[vn][0,-1,:,:].values
            vmin = -2
            vmax = 2
            cmap='hsv_r'
        elif vn == 'zeta':
            v = ds[vn][0,:,:].values
            h = ds.h.values
            mr = ds.mask_rho.values
            v[mr==0] = np.nan
            h[mr==0] = np.nan
            v = v + h
            vn = 'depth'
            vmin = 2
            vmax = 4
            cmap='RdYlGn'
        else:
            v = ds[vn][0, -1,:,:].values
        ax = fig.add_subplot(1, len(vn_list), ii)
        ax.set_xticks([])
        ax.set_yticks([])
        cs = ax.pcolormesh(px, py, v, cmap=cmap, vmin=vmin, vmax=vmax)
        #pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'], his_num=True)
        vmax, vjmax, vimax, vmin, vjmin, vimin = pfun.maxmin(v)
        ax.plot(x[vjmax,vimax], y[vjmax,vimax],'*y', mec='k', markersize=15)
        ax.plot(x[vjmin,vimin], y[vjmin,vimin],'oy', mec='k', markersize=10)
        ax.set_title(('%s ((*)max=%0.1f, (o)min=%0.1f)' % (vn, vmax, vmin)))
        ii += 1

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_basic(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    fig = plt.figure()
    # PLOT CODE
    vn_list = ['salt', 'temp']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], vlims_fac=pinfo.range_dict[vn])
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            ax.set_yticklabels([])
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_Fb(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    fig = plt.figure()
    # PLOT CODE
    vn_list = ['salt', 'Fb']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        if vn == 'salt':
            cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                    cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                    vlims_fac=pinfo.range_dict[vn])
            ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]),
                    fontsize=1.2*fs)
        elif vn == 'Fb':
            # plot vertically integrateed buoyancy flux
            # calculate potential density
            import seawater as sw
            rho = sw.dens0(ds.salt.values.squeeze(), ds.temp.values.squeeze())
            # calculate vertically-integrated buoyancy flux
            Fb = -9.8 * np.sum(ds.AKs[0,1:-1,:,:].squeeze() * np.diff(rho, axis=0), axis=0).values
            Fb[ds.mask_rho.values==0] = np.nan
            plon, plat = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
            cs = ax.pcolormesh(plon,plat, Fb, vmin=0, vmax=.5, cmap='YlOrRd')
            ax.set_title('Vertically Integrated Fb [W/m2]')
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
            #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        elif ii == 2:
            ax.set_yticklabels([])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_fancy(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    fig = plt.figure()
    # PLOT CODE
    vn_list = ['salt', 'temp']
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        if vn == 'salt':
            cmap = 'jet'
            vlims_fac = .5
        elif vn == 'temp':
            cmap = 'RdYlBu_r'
            vlims_fac = 1
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
                cmap=cmap, fac=pinfo.fac_dict[vn], vlims_fac=vlims_fac)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            ax.set_yticklabels([])
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_dive_vort(in_dict):
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 10
    pfun.start_plot(fs=fs, figsize=(int(hgt*2.5/AR),int(hgt)))
    fig = plt.figure()
    
    # create fields
    u = ds.u[0,-1,:,:].values
    v = ds.v[0,-1,:,:].values
    dx = 1/ds.pm.values
    dy = 1/ds.pn.values
    # dive is on the trimmed rho grid
    dive = np.diff(u[1:-1,:], axis=1)/dx[1:-1,1:-1] + np.diff(v[:,1:-1],axis=0)/dy[1:-1,1:-1]
    # vort is on the psi grid (plot with lon_rho, lat_rho)
    vort = np.diff(v,axis=1)/dx[1:,1:] - np.diff(u,axis=0)/dy[1:,1:]
    
    # set color limits
    vv = 2*np.nanstd(vort)
    
    # PLOT CODE
    if in_dict['auto_vlims']:
        pinfo.vlims_dict['vort'] = (-vv, vv)
        pinfo.vlims_dict['dive'] = (-vv, vv)
        
    vmin = pinfo.vlims_dict['vort'][0]
    vmax = pinfo.vlims_dict['vort'][1]
    
    for ii in [1,2]:
        ax = fig.add_subplot(1, 2, ii)
        cmap = 'RdYlBu_r'
        if ii == 1:
            plon, plat = pfun.get_plon_plat(ds.lon_rho[1:-1,1:-1].values, ds.lat_rho[1:-1,1:-1].values)
            cs = plt.pcolormesh(plon, plat, dive, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Divergence $[s^{-1}]$', fontsize=1.2*fs)
        elif ii == 2:
            cs = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, vort, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Vorticity $[s^{-1}]$', fontsize=1.2*fs)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            pass
            #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_dive_vort2(in_dict):
    # same as dive_vort but focused on a specific region
    
    # JdF:
    aa = [-125, -122.3, 47.8, 48.8]
    
    # START
    ds = xr.open_dataset(in_dict['fn'])
    # find aspect ratio of the map
    # aa = pfun.get_aa(ds)
    # AR is the aspect ratio of the map: Vertical/Horizontal
    AR = (aa[3] - aa[2]) / (np.sin(np.pi*aa[2]/180)*(aa[1] - aa[0]))
    fs = 14
    hgt = 6
    pfun.start_plot(fs=fs, figsize=(10,10))
    fig = plt.figure()
    
    # create fields
    u = ds.u[0,-1,:,:].values
    v = ds.v[0,-1,:,:].values
    dx = 1/ds.pm.values
    dy = 1/ds.pn.values
    # dive is on the trimmed rho grid
    dive = np.diff(u[1:-1,:], axis=1)/dx[1:-1,1:-1] + np.diff(v[:,1:-1],axis=0)/dy[1:-1,1:-1]
    # vort is on the psi grid (plot with lon_rho, lat_rho)
    vort = np.diff(v,axis=1)/dx[1:,1:] - np.diff(u,axis=0)/dy[1:,1:]
    
    # set color limits
    vv = 4*np.nanstd(vort)
    
    # PLOT CODE
    if in_dict['auto_vlims']:
        pinfo.vlims_dict['vort'] = (-vv, vv)
        pinfo.vlims_dict['dive'] = (-vv, vv)
        
    vmin = pinfo.vlims_dict['vort'][0]
    vmax = pinfo.vlims_dict['vort'][1]
    
    for ii in [1,2]:
        ax = fig.add_subplot(2, 1, ii)
        cmap = 'RdYlBu_r'
        if ii == 1:
            plon, plat = pfun.get_plon_plat(ds.lon_rho[1:-1,1:-1].values, ds.lat_rho[1:-1,1:-1].values)
            cs = plt.pcolormesh(plon, plat, dive, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Divergence $[s^{-1}]$', fontsize=1.2*fs)
        elif ii == 2:
            cs = plt.pcolormesh(ds.lon_rho.values, ds.lat_rho.values, vort, cmap=cmap, vmin = vmin, vmax = vmax)
            ax.set_title('Surface Vorticity $[s^{-1}]$', fontsize=1.2*fs)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(aa)
        pfun.dar(ax)
        ax.set_ylabel('Latitude')
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'])
            #pfun.add_windstress_flower(ax, ds)
            #pfun.add_bathy_contours(ax, ds, txt=True)
        elif ii == 2:
            ax.set_xlabel('Longitude')
            #pfun.add_velocity_vectors(ax, ds, in_dict['fn'])
        ii += 1
    #fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_ri(in_dict):
    """
    Simplified Richardson number
    """
    # START
    fs = 10
    pfun.start_plot(fs=fs, figsize=(20,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PLOT CODE
    xrho = ds['lon_rho'][0,:].values
    yrho = ds['lat_rho'][:,0].values

    # define box
    aa = [-123.25, -122.1, 47, 48.75]
    ix0 = zfun.find_nearest_ind(xrho, aa[0])
    ix1 = zfun.find_nearest_ind(xrho, aa[1])
    iy0 = zfun.find_nearest_ind(yrho, aa[2])
    iy1 = zfun.find_nearest_ind(yrho, aa[3])

    h = ds.h[iy0:iy1, ix0:ix1].values
    rho_bot = ds.rho[0, 0, iy0:iy1, ix0:ix1].values
    rho_top = ds.rho[0, -1, iy0:iy1, ix0:ix1].values
    drho = rho_bot - rho_top
    u = ds.ubar[0, iy0:iy1, ix0-1:ix1].values
    v = ds.vbar[0, iy0-1:iy1, ix0:ix1].values
    u[np.isnan(u)] = 0
    v[np.isnan(v)] = 0
    uu = (u[:, 1:] + u[:, :-1])/2
    vv = (v[1:, :] + v[:-1, :])/2
    spd2 = uu**2 + vv**2
    spd2[np.isnan(drho)] = np.nan
    spd2[spd2 < .001] = .001 # avoid divide by zero errors

    # approximate Richardson number
    rho0 = ds.rho0.values
    g = 9.8
    Ri = g * drho * h / (rho0 * spd2)

    # psi_grid coordinates
    x, y = np.meshgrid(ds.lon_u.values[0,ix0-1:ix1], ds.lat_v.values[iy0-1:iy1,0])

    # PLOTTING
    plt.close('all')
    pfun.start_plot(fs=10, figsize=(18,10))
    fig = plt.figure()

    xt = [-123.2, -122.2]
    yt = [47, 47.5, 48, 48.5]

    ax = fig.add_subplot(131)
    cs = ax.pcolormesh(x, y, drho, vmin=0, vmax=5, cmap=cm.dense)
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    ax.set_title(r'$\Delta\rho\ [kg\ m^{-3}]$')
    ax.set_xticks(xt)
    ax.set_yticks(yt)

    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(x, y, np.sqrt(spd2), vmin=0, vmax=2, cmap=cm.speed)
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    ax.set_title(r'Speed $[m\ s^{-1}]$')
    ax.set_xticks(xt)
    ax.set_yticks(yt)
    ax.set_yticklabels([])

    ax = fig.add_subplot(133)
    cs = ax.pcolormesh(x, y, 4*Ri, vmin=0, vmax = 2, cmap='RdYlBu')
    fig.colorbar(cs, ax=ax)
    pfun.dar(ax)
    pfun.add_coast(ax)
    ax.axis(aa)
    ax.set_title(r'$4 x Ri$')
    ax.set_xticks(xt)
    ax.set_yticks(yt)
    ax.set_yticklabels([])
        
    fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_Chl_DO(in_dict):
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn_list = ['phytoplankton', 'oxygen']
    fs = 14
    ii = 1
    for vn in vn_list:
        if vn == 'phytoplankton':
            slev = -1
            stext = 'Surface'
        elif vn == 'oxygen':
            slev = 0
            stext = 'Bottom'
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
                cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
                vlims_fac=pinfo.range_dict[vn], do_mask_edges=True)
        fig.colorbar(cs)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_title('%s %s %s' % (stext, pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
        ax.set_xlabel('Longitude')
        pfun.add_bathy_contours(ax, ds, txt=True)
        if ii == 1:
            ax.set_ylabel('Latitude')
            pfun.add_info(ax, in_dict['fn'])
            pfun.add_windstress_flower(ax, ds)
        ii += 1
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_DO_WA_shelf(in_dict):
    # Focus on bottom DO on the WA shelf
    aa = [-126.1, -123.7, 45.8, 48.8]
    xtl = [-126, -125, -124]
    ytl = [46, 47, 48]
    
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(7,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'oxygen'
    slev = 0
    stext = 'Bottom'
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = ()
    ax = fig.add_subplot(111)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict, slev=slev,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn],
            vlims_fac=pinfo.range_dict[vn], do_mask_edges=True)
    fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    ax.set_title('%s %s %s' % (stext, pinfo.tstr_dict[vn],pinfo.units_dict[vn]), fontsize=1.2*fs)
    ax.set_xlabel('Longitude')
    pfun.add_bathy_contours(ax, ds, txt=False)
    ax.set_ylabel('Latitude')
    ax.set_xticks(xtl)
    ax.set_yticks(ytl)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    
    pfun.add_windstress_flower(ax, ds, t_scl=0.5, t_leglen=0.1, center=(.85,.65), fs=12)
    # ADD MEAN WINDSTRESS VECTOR
    # t_scl: scale windstress vector (smaller to get longer arrows)
    # t_leglen: # Pa for wind stress vector legend
        
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_ths(in_dict):
    # Plot property-property plots, like theta vs. s
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(10,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    # make a potential density field
    import seawater as sw
    s0 = 25; s1 = 35
    th0 = 0; th1 = 20
    SS, TH = np.meshgrid(np.linspace(s0, s1, 50), np.linspace(th0, th1, 50))
    SIG = sw.dens0(SS, TH) - 1000
    S = zrfun.get_basic_info(in_dict['fn'], only_S=True)
    h = ds['h'].values
    z = zrfun.get_z(h, 0*h, S, only_rho=True)
    s = ds['salt'].values.squeeze()
    th = ds['temp'].values.squeeze()
    ax = fig.add_subplot(111)
    ax.set_xlabel('Salinity')
    ax.set_ylabel('Theta (deg C)')
    ax.contour(SS, TH, SIG, 20)
    nsub = 500
    alpha = .1
    mask = z > -10
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.r', alpha=alpha)
    mask = (z < -10) & (z > -200)
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.g', alpha=alpha)
    mask = z < -200
    ax.plot(s[mask][::nsub], th[mask][::nsub], '.b', alpha=alpha)
    ax.set_xlim(s0, s1)
    ax.set_ylim(th0, th1)
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_debug(in_dict):
    # Focused on debugging
    vn_list = ['u', 'v', 'zeta']
    do_wetdry = False
    
    # START
    fs = 10
    pfun.start_plot(fs=fs, figsize=(8*len(vn_list),10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    ii = 1
    for vn in vn_list:
        if 'lon_rho' in ds[vn].coords:
            tag = 'rho'
        if 'lon_u' in ds[vn].coords:
            tag = 'u'
        if 'lon_v' in ds[vn].coords:
            tag = 'v'
        x = ds['lon_'+tag].values
        y = ds['lat_'+tag].values
        px, py = pfun.get_plon_plat(x,y)
        if vn in ['u', 'v']:
            v = ds[vn][0,-1,:,:].values
            vmin = -2
            vmax = 2
            cmap='hsv_r'
        elif vn == 'zeta':
            v = ds[vn][0,:,:].values
            h = ds.h.values
            mr = ds.mask_rho.values
            v[mr==0] = np.nan
            h[mr==0] = np.nan
            v = v + h
            vn = 'depth'
            vmin = 2
            vmax = 4
            cmap='RdYlGn'
        else:
            v = ds[vn][0, -1,:,:].values
        ax = fig.add_subplot(1, len(vn_list), ii)
        ax.set_xticks([])
        ax.set_yticks([])
        cs = ax.pcolormesh(px, py, v, cmap=cmap, vmin=vmin, vmax=vmax)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'], his_num=True)
        vmax, vjmax, vimax, vmin, vjmin, vimin = pfun.maxmin(v)
        ax.plot(x[vjmax,vimax], y[vjmax,vimax],'*y', mec='k', markersize=15)
        ax.plot(x[vjmin,vimin], y[vjmin,vimin],'oy', mec='k', markersize=10)
        ax.set_title(('%s ((*)max=%0.1f, (o)min=%0.1f)' % (vn, vmax, vmin)))
        ii += 1

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_layer(in_dict):
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn_list = ['oxygen', 'temp']
    z_level = -250
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, z_level)
        v_scaled = pinfo.fac_dict[vn]*laym
        vlims = pinfo.vlims_dict[vn]
        if len(vlims) == 0:
            vlims = pfun.auto_lims(v_scaled)
            pinfo.vlims_dict[vn] = vlims
        cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], v_scaled[1:-1,1:-1],
                           vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
        cb = fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        ax.set_title('%s %s on Z = %d (m)' % (pinfo.tstr_dict[vn], pinfo.units_dict[vn], z_level))
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'])
            ax.set_ylabel('Latitude')
            pfun.add_windstress_flower(ax, ds)
        if ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], zlev=z_level)
        ii += 1
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_bpress(in_dict):
    """
    Specialized plot related to bottom pressure anomalies.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(14,10))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    sta_dict = {
        'CE01':(-124.095, 44.6598), # Oregon Inshore (25 m)
        'CE02':(-124.304, 44.6393), # Oregon Shelf (80 m)
        'CE04':(-124.956, 44.3811), # Oregon Offshore (588 m)
        'PN01A':(-125.3983, 44.5096), # Slope Base (2905 m)
        }
    vn_list = ['salt', 'temp']
    z_level = -300
    zfull = pfun.get_zfull(ds, in_dict['fn'], 'rho')
    ii = 1
    for vn in vn_list:
        if in_dict['auto_vlims']:
            pinfo.vlims_dict[vn] = ()
        ax = fig.add_subplot(1, len(vn_list), ii)
        laym = pfun.get_laym(ds, zfull, ds['mask_rho'][:], vn, z_level)
        v_scaled = pinfo.fac_dict[vn]*laym
        vlims = pinfo.vlims_dict[vn]
        if len(vlims) == 0:
            if vn == 'salt':
                vlims = pfun.auto_lims(v_scaled, vlims_fac=0.3)
            elif vn == 'temp':
                vlims = pfun.auto_lims(v_scaled, vlims_fac=2)
            else:
                vlims = pfun.auto_lims(v_scaled)
            pinfo.vlims_dict[vn] = vlims
        cs = ax.pcolormesh(ds['lon_psi'][:], ds['lat_psi'][:], v_scaled[1:-1,1:-1],
                           vmin=vlims[0], vmax=vlims[1], cmap=pinfo.cmap_dict[vn])
        cb = fig.colorbar(cs)
        pfun.add_bathy_contours(ax, ds, txt=True)
        pfun.add_coast(ax)
        ax.axis(pfun.get_aa(ds))
        pfun.dar(ax)
        ax.set_xlabel('Longitude')
        ax.set_title('%s %s on Z = %d (m)' % (pinfo.tstr_dict[vn], pinfo.units_dict[vn], z_level))
        if ii == 1:
            pfun.add_info(ax, in_dict['fn'])
            ax.set_ylabel('Latitude')
            pfun.add_windstress_flower(ax, ds)
        if ii == 2:
            pfun.add_velocity_vectors(ax, ds, in_dict['fn'], zlev=z_level, v_scl=5, v_leglen=0.3)
            for sta in sta_dict.keys():
                ax.plot(sta_dict[sta][0], sta_dict[sta][1], 'o', mfc='y', mec='k', ms=10)
        ii += 1
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_sect(in_dict):
    """
    This plots a map and a section (distance, z), and makes sure
    that the color limits are identical.  If the color limits are
    set automatically then the section is the preferred field for
    setting the limits.
    
    I think this works best with -avl False (the default).
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(20,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    # PLOT CODE
    vn = 'phytoplankton'
    # GET DATA
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # CREATE THE SECTION
    # create track by hand
    if False:
        lon = G['lon_rho']
        lat = G['lat_rho']
        zdeep = -3500
        x = np.linspace(lon.min(), lon.max(), 500)
        y = 47 * np.ones(x.shape)
    # or read in a section (or list of sections)
    else:
        tracks_path = Ldir['data'] / 'section_lines'
        tracks = ['Line_jdf_v0.p', 'Line_ps_main_v0.p']
        zdeep = -300
        xx = np.array([])
        yy = np.array([])
        for track in tracks:
            track_fn = tracks_path / track
            # get the track to interpolate onto
            pdict = pickle.load(open(track_fn, 'rb'))
            xx = np.concatenate((xx,pdict['lon_poly']))
            yy = np.concatenate((yy,pdict['lat_poly']))
        for ii in range(len(xx)-1):
            x0 = xx[ii]
            x1 = xx[ii+1]
            y0 = yy[ii]
            y1 = yy[ii+1]
            nn = 20
            if ii == 0:
                x = np.linspace(x0, x1, nn)
                y = np.linspace(y0,y1, nn)
            else:
                x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
                y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)
    
    # COLOR
    # scaled section data
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # now we use the scaled section as the preferred field for setting the
    # color limits of both figures in the case -avl True
    if in_dict['auto_vlims']:
        pinfo.vlims_dict[vn] = pfun.auto_lims(sf)
    
    # PLOTTING
    # map with section line
    ax = fig.add_subplot(1, 3, 1)
    cs = pfun.add_map_field(ax, ds, vn, pinfo.vlims_dict,
            cmap=pinfo.cmap_dict[vn], fac=pinfo.fac_dict[vn], do_mask_edges=True)
    # fig.colorbar(cs, ax=ax) # It is identical to that of the section
    pfun.add_coast(ax)
    aaf = [-125.5, -122.1, 46.8, 50.3] # focus domain
    ax.axis(aaf)
    pfun.dar(ax)
    pfun.add_info(ax, in_dict['fn'], loc='upper_right')
    ax.set_title('Surface %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    # add section track
    ax.plot(x, y, '-r', linewidth=2)
    ax.plot(x[idist0], y[idist0], 'or', markersize=5, markerfacecolor='w',
        markeredgecolor='r', markeredgewidth=2)
    ax.set_xticks([-125, -124, -123])
    ax.set_yticks([47, 48, 49, 50])
    # section
    ax = fig.add_subplot(1, 3, (2, 3))
    ax.plot(dist, v2['zbot'], '-k', linewidth=2)
    ax.plot(dist, v2['zeta'], '-b', linewidth=1)
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 5)
    # plot section
    svlims = pinfo.vlims_dict[vn]
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=svlims[0], vmax=svlims[1], cmap=pinfo.cmap_dict[vn])
    fig.colorbar(cs, ax=ax)
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Z (m)')
    ax.set_title('Section %s %s' % (pinfo.tstr_dict[vn],pinfo.units_dict[vn]))
    fig.tight_layout()
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_sect_soundspeed(in_dict):
    """
    Soundspeed section plot
    """
    import gsw

    ds = xr.open_dataset(in_dict['fn'])
    # create track by hand
    x = np.linspace(-124.85,-124.2, 100) # shelf only
    #x = np.linspace(-126,-124.2, 100) # shows SOFAR channel
    y = 47 * np.ones(x.shape)
    v2, v3, dist, idist0 = pfun.get_section(ds, 'salt', x, y, in_dict)
    s = v3['sectvarf']
    v2, v3, dist, idist0 = pfun.get_section(ds, 'temp', x, y, in_dict)
    th = v3['sectvarf']

    X = v3['distf']
    Z = v3['zrf']
    # adjust Z so surface is at 0
    Z = Z - Z[-1,:]

    p = gsw.p_from_z(Z, 47)
    SA = gsw.SA_from_SP(s, p, -125, 47)
    CT = gsw.CT_from_pt(SA, th)
    spd = gsw.sound_speed(SA, CT, p)

    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(16,9))
    fig, axes = plt.subplots(nrows=3, ncols=2)

    ax = axes[0,0]
    cs = ax.pcolormesh(X, Z, SA, cmap='jet')
    fig.colorbar(cs, ax=ax)
    ax.text(.95, .05, 'Absolute Salinity', transform=ax.transAxes, ha='right')

    ax = axes[1,0]
    cs = ax.pcolormesh(X, Z, CT, cmap='jet')
    fig.colorbar(cs, ax=ax)
    ax.text(.95, .05, 'Conservative Temperature', transform=ax.transAxes, ha='right')

    ax = axes[2,0]
    cs = ax.pcolormesh(X, Z, spd, cmap='jet')
    fig.colorbar(cs, ax=ax)
    ax.text(.95, .05, 'Soundspeed [m/s]', transform=ax.transAxes, ha='right')

    ax = axes[0,1]
    ax.plot(SA,Z, alpha=.2)
    ax.text(.05, .05, 'Absolute Salinity', transform=ax.transAxes, ha='left')

    ax = axes[1,1]
    ax.plot(CT,Z, alpha=.2)
    ax.text(.95, .05, 'Conservative Temperature', transform=ax.transAxes, ha='right')

    ax = axes[2,1]
    ax.plot(spd,Z, alpha=.2)
    ax.text(.95, .05, 'Soundspeed [m/s]', transform=ax.transAxes, ha='right')

    fig.suptitle(str(in_dict['fn']))

    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
    

def P_splash(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  Eventually I could automate making this new every day.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from PyCO2SYS import CO2SYS
    import seawater as sw
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from warnings import filterwarnings
    filterwarnings('ignore') # skip a warning from PyCO2SYS

    Ldir = Lfun.Lstart()

    do_topo = True

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x = ds.lon_psi.values
    y = ds.lat_psi.values
    th = ds['temp'][0,-1,1:-1,1:-1].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 3
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    def get_arag(ds, fn, aa, nlev):
        G = zrfun.get_basic_info(fn, only_G=True)
        # find indices that encompass region aa
        i0 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[0]) - 1
        i1 = zfun.find_nearest_ind(G['lon_rho'][0,:], aa[1]) + 2
        j0 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[2]) - 1
        j1 = zfun.find_nearest_ind(G['lat_rho'][:,0], aa[3]) + 2
        px = G['lon_psi'][j0:j1-1, i0:i1-1]
        py = G['lat_psi'][j0:j1-1, i0:i1-1]
        lat = G['lat_rho'][j0:j1,i0:i1] # used in sw.pres
        # first extract needed fields and save in v_dict
        v_dict = {}
        vn_in_list = ['temp', 'salt' , 'rho', 'alkalinity', 'TIC']
        for cvn in vn_in_list:
            L = ds[cvn][0,nlev,j0:j1,i0:i1].values
            v_dict[cvn] = L
        # ------------- the CO2SYS steps -------------------------
        # create pressure
        Ld = G['h'][j0:j1,i0:i1]
        Lpres = sw.pres(Ld, lat)
        # get in situ temperature from potential temperature
        Ltemp = sw.ptmp(v_dict['salt'], v_dict['temp'], 0, Lpres)
        # convert from umol/L to umol/kg using in situ dentity
        Lalkalinity = 1000 * v_dict['alkalinity'] / (v_dict['rho'] + 1000)
        Lalkalinity[Lalkalinity < 100] = np.nan
        LTIC = 1000 * v_dict['TIC'] / (v_dict['rho'] + 1000)
        LTIC[LTIC < 100] = np.nan
        CO2dict = CO2SYS(Lalkalinity, LTIC, 1, 2, v_dict['salt'], Ltemp, Ltemp,
            Lpres, Lpres, 50, 2, 1, 10, 1, NH3=0.0, H2S=0.0)
        # PH = CO2dict['pHout']
        # PH = zfun.fillit(PH.reshape((v_dict['salt'].shape)))
        ARAG = CO2dict['OmegaARout']
        ARAG = ARAG.reshape((v_dict['salt'].shape))
        ARAG = ARAG[1:-1, 1:-1]
        return px, py, ARAG

    # LARGE MAP
    ax = fig.add_subplot(121)
    cmap = 'RdYlBu_r'
    cs = ax.pcolormesh(x,y,th, cmap=cmap, vmin=11, vmax=20)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmap = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmap, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.2,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.98,.95,'Surface water\nTemperature $[^{\circ}C]$\n'+tstr,
         ha='right', va='top', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Willapa and Grays Harbor
    aa = [-124.6, -123.65, 46, 47.2]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='g', alpha=1, linewidth=2, inset=0)

    # SMALL. MAP
    ax = fig.add_subplot(122)
    px, py, ARAG = get_arag(ds, fn, aa, nlev)
    cs = ax.pcolormesh(px,py,ARAG, cmap='coolwarm_r', vmin=0, vmax=3)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="4%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-124.5, -124])
    ax.set_yticks([46, 47])
    ax.set_xlabel('Longitude')
    ax.text(.98,.99,'Bottom water\nAragonite\nSaturation\nState',
         ha='right', va='top', weight='bold', transform=ax.transAxes)

    fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_splash2(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  This one is focused on the Salish Sea.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(12,9))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    Ldir = Lfun.Lstart()

    do_topo = True

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x = ds.lon_psi.values
    y = ds.lat_psi.values
    th = ds['temp'][0,-1,1:-1,1:-1].values
    ox = pinfo.fac_dict['oxygen'] * ds['oxygen'][0,0,1:-1,1:-1].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 3
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    # LARGE MAP
    ax = fig.add_subplot(121)
    cmap = 'RdYlBu_r'
    cs = ax.pcolormesh(x,y,th, cmap=cmap, vmin=11, vmax=20)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmap = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmap, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.5,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.03,.45,'Surface water\nTemperature $[^{\circ}C]$\n'+tstr,
         ha='left', va='bottom', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Salish Sea
    aa = [-125.3, -122.1, 46.8, 50.3]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='g', alpha=1, linewidth=2, inset=0)
    
    fs2 = fs*.9
    fs3 = fs*.8
    
    ax.text(-123.072,46.7866,'Washington', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-122.996,44.5788,'Oregon', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    
    ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
        style='italic',ha='center',va='center',rotation=-10)
    

    # SMALL MAP
    from cmocean import cm
    
    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x,y,ox, cmap=cm.oxy, vmin=0, vmax=10)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="30%", loc='upper right', borderpad=2)
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-125, -124, -123])
    ax.set_yticks([47, 48, 49, 50])
    ax.set_xlabel('Longitude')
    ax.text(.84,.95,'Salish Sea\n\nBottom Oxygen\n$[mg\ L^{-1}]$',
         ha='right', va='top', weight='bold', transform=ax.transAxes)
         
    # add labels
    ax.text(-122.8,49.335,'Fraser\nRiver',size=fs2,
        style='italic',ha='center',va='center',rotation=0)
    ax.text(-123.7,49.2528,'Strait of Georgia',size=fs2,
        style='italic',ha='center',va='center',rotation=-30)
    ax.text(-123.5,48.28,'Strait of Juan de Fuca',size=fs2,
        style='italic',ha='center',va='center',rotation=0,
        color='w')
    ax.text(-123.3,47.6143,'Puget\nSound',size=fs2,
        style='italic',ha='center',va='center',rotation=+55)
    ax.text(-122.3,48.48,'Skagit\nRiver',size=fs3,
        style='italic',ha='center',va='center',
        bbox=dict(facecolor='w', edgecolor='None',alpha=.5))
    ax.text(-123.173,48.44,'Haro\nStrait',size=fs3,
        style='italic',ha='center',va='center',
        color='w')

    fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_splash3(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  This one is focused on the Puget Sound.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(15,12))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    Ldir = Lfun.Lstart()

    do_topo = True

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x,y = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
    th = ds['temp'][0,-1,:,:].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 3
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    # LARGE MAP
    ax = fig.add_subplot(121)
    cs = ax.pcolormesh(x,y,th, cmap='RdYlBu_r', vmin=11, vmax=20)
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmap = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap='gist_earth', shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.5,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.03,.45,'Surface water\nTemperature $[^{\circ}C]$\n'+tstr,
         ha='left', va='bottom', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Puget Sound and the San Juans
    aa = [-123.4, -122, 47, 48.8]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='m', alpha=1, linewidth=2, inset=0)
    
    fs2 = fs*.9
    fs3 = fs*.8
    
    ax.text(-123.072,46.7866,'Washington', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-122.996,44.5788,'Oregon', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    
    ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
        style='italic',ha='center',va='center',rotation=-10)
    

    # SMALL MAP
    
    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x,y,th, cmap='RdYlBu_r', vmin=11, vmax=20)
    # Inset colorbar
    # cbaxes = inset_axes(ax, width="5%", height="30%", loc='upper right', borderpad=2)
    # fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-123, -122.5, -122])
    ax.set_yticks([47, 48])
    ax.set_xlabel('Longitude')
    ax.text(.03,.5,'Puget Sound &\nSan Juans',
         ha='left', va='center', weight='bold', transform=ax.transAxes)
             

    #fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
        
def P_splash4(in_dict):
    """
    This makes a fancy plot suitable for the landing page of the LiveOcean
    website.  This one is focused on the Puget Sound Salinity.
    """
    # START
    fs = 14
    pfun.start_plot(fs=fs, figsize=(15,12))
    fig = plt.figure()
    ds = xr.open_dataset(in_dict['fn'])
    
    # PREPARING FIELDS
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes

    Ldir = Lfun.Lstart()

    do_topo = True
    cmap = 'Spectral_r'
    vlims = (25, 33) # full map
    vlims2 = (25, 33) # PS map

    # model output
    fn = in_dict['fn']
    T = zrfun.get_basic_info(fn, only_T=True)
    x,y = pfun.get_plon_plat(ds.lon_rho.values, ds.lat_rho.values)
    salt = ds['salt'][0,-1,:,:].values

    if do_topo:
        # topography
        tfn = (Ldir['data'] / 'topo' / 'srtm15' / 'topo15.nc')
        tds = xr.open_dataset(tfn)
        step = 1
        tx = tds['lon'][::step].values
        ty = tds['lat'][::step].values
        tz = tds['z'][::step,::step].values
        tz[tz<0] = np.nan

    # LARGE MAP
    ax = fig.add_subplot(121)
    cs = ax.pcolormesh(x,y,salt, cmap=cmap, vmin=vlims[0], vmax=vlims[1])
    # Inset colorbar
    cbaxes = inset_axes(ax, width="5%", height="40%", loc='lower left')
    fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmapt = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmapt, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis([-130, -122, 42, 52])
    ax.set_xticks([-129, -127, -125, -123])
    ax.set_yticks([42, 44, 46, 48, 50, 52])
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    tstr = T['dt'].strftime(Lfun.ds_fmt)
    ax.text(.98,.99,'LiveOcean', size=fs*1.5,
         ha='right', va='top', weight='bold', transform=ax.transAxes)
    ax.text(.03,.45,'Surface Salinity $[g\ kg^{-1}]$\n'+tstr,
         ha='left', va='bottom', weight='bold', transform=ax.transAxes,
         bbox=dict(facecolor='w', edgecolor='None',alpha=.5))

    # box for Puget Sound and the San Juans
    aa = [-123.4, -122, 47, 48.8]
    nlev = 0
    # draw box on the large map
    pfun.draw_box(ax, aa, linestyle='-', color='m', alpha=1, linewidth=2, inset=0)
    
    fs2 = fs*.9
    fs3 = fs*.8
    
    ax.text(-123.072,46.7866,'Washington', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-122.996,44.5788,'Oregon', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    
    ah = ax.text(-125.3,49.4768,'Vancouver\nIsland', size=fs2,
        style='italic',ha='center',va='center',rotation=-45)
    ax.text(-126.3,50.2,'Johnstone\nStrait', size=.7*fs2,
        style='italic',ha='center',va='center',rotation=-10)
    

    # SMALL MAP
    
    ax = fig.add_subplot(122)
    cs = ax.pcolormesh(x,y,salt, cmap=cmap, vmin=vlims2[0], vmax=vlims2[1])
    # Inset colorbar
    # cbaxes = inset_axes(ax, width="5%", height="30%", loc='upper right', borderpad=2)
    # fig.colorbar(cs, cax=cbaxes, orientation='vertical')
    if do_topo:
        cmapt = 'gist_earth'
        cs = ax.pcolormesh(tx,ty,tz, cmap=cmapt, shading='nearest', vmin=-1000, vmax=2000)
    pfun.add_coast(ax)
    pfun.dar(ax)
    ax.axis(aa)
    ax.set_xticks([-123, -122.5, -122])
    ax.set_yticks([47, 48])
    ax.set_xlabel('Longitude')
    ax.text(.03,.5,'Puget Sound &\nSan Juans',
         ha='left', va='center', weight='bold', transform=ax.transAxes)
             

    #fig.tight_layout()
    
    # FINISH
    ds.close()
    pfun.end_plot()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_salt(in_dict):
    # Plot salinity maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.

    vn = 'salt'
    vlims = (28.5, 33) # full map
    vlims2 = (22, 32) # PS map
    vlims3 = (29, 32) # PS section
    cmap = 'Spectral_r'

    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] / 'section_lines'
    tracks = ['Line_ps_main_v0.p']
    zdeep = -300
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path / track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))
    fs = 18 # fontsize

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'].values
    lat = ds['lat_psi'].values
    v =ds[vn][0, -1, 1:-1, 1:-1].values
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nSalinity\n'
        + datetime.strftime(T['dt'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')
    
    ax.text(.99,.97,'S range\n'+ str(vlims), transform=ax.transAxes,
        va='top', ha='right', c='orange', size=.6*fs, weight='bold')

    # PS map
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims2[0], vmax=vlims2[1],
        cmap=cmap)
    #fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    sect_color = 'violet'
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle='--', color='k', linewidth=2)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')
    ax.text(.93,.97,'S range\n'+ str(vlims2), transform=ax.transAxes,
        va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    

    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, v2['zeta']+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    ax.text(.99,.4,'S range\n'+ str(vlims3), transform=ax.transAxes,
        va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # Tides
    alpha = .4
    ax = fig.add_subplot(436)
    ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
        lw=3, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(0,1.5)
    ax.set_axis_off()

    # Wind
    alpha=.5
    ax = fig.add_subplot(439)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
        color='g', alpha=alpha, fontsize=fs)
    ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
        color='b', alpha=alpha, fontsize=fs)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.15, .25)
    ax.set_axis_off()

    # Rivers
    alpha = .6
    cr = fdf['Columbia R. Flow (1000 m3/s)'].values
    fr = fdf['Fraser R. Flow (1000 m3/s)'].values
    sr = fdf['Skagit R. Flow (1000 m3/s)'].values
    this_yd = fdf.loc[TM, 'yearday']
    ax = fig.add_subplot(4,3,12)
    ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
    ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
    ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)
    # time markers
    ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(.9, .85, 'Columbia River', transform=ax.transAxes,
        color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .70, 'Fraser River', transform=ax.transAxes,
        color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .55, 'Skagit River', transform=ax.transAxes,
        color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-5,20)
    ax.set_axis_off()

    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = 'r'
    if tm.month in [4, 5, 6]:
        clist[1] = 'r'
    if tm.month in [7, 8, 9]:
        clist[2] = 'r'
    if tm.month in [10, 11, 12]:
        clist[3] = 'r'
    ax.text(0, 0, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=fs, horizontalalignment='left', style='italic')
    ax.text(.4, 0, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(.68, 0, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(1, 0, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=fs, horizontalalignment='right', style='italic')

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_oxygen(in_dict):
    # Plot bottom oxygen maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.

    vn = 'oxygen'
    vlims = (0, 10) # full map
    vlims2 = (0, 10) # PS map
    vlims3 = (0, 10) # PS section
    from cmocean import cm
    cmap = cm.oxy

    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] / 'section_lines'
    tracks = ['Line_HC_thalweg_long.p']
    zdeep = -250
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path / track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))
    fs = 18 # fontsize

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'].values
    lat = ds['lat_psi'].values
    v =ds[vn][0, 0, 1:-1, 1:-1].values
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nBottom Oxygen\n'
        + datetime.strftime(T['dt'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cbaxes = inset_axes(ax, width="40%", height="4%", loc='upper right', borderpad=2) 
    cb = fig.colorbar(cs, cax=cbaxes, orientation='horizontal')
    cb.ax.tick_params(labelsize=.85*fs)
    ax.text(1, .85, r'$[mg\ L^{-1}]$', transform=ax.transAxes, fontsize=fs, ha='right')
    # ax.text(.99,.97,'S range\n'+ str(vlims), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')

    # PS map
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims2[0], vmax=vlims2[1],
        cmap=cmap)
    #fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    sect_color = 'violet'
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle='--', color='k', linewidth=2)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')
    # ax.text(.93,.97,'S range\n'+ str(vlims2), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    

    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, v2['zeta']+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    # ax.text(.99,.4,'S range\n'+ str(vlims3), transform=ax.transAxes,
    #     va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nHood Canal', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # Tides
    alpha = .4
    ax = fig.add_subplot(436)
    ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
        lw=3, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(0,1.5)
    ax.set_axis_off()

    # Wind
    alpha=.5
    ax = fig.add_subplot(439)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
        color='g', alpha=alpha, fontsize=fs)
    ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
        color='b', alpha=alpha, fontsize=fs)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.15, .25)
    ax.set_axis_off()

    # Rivers
    alpha = .6
    cr = fdf['Columbia R. Flow (1000 m3/s)'].values
    fr = fdf['Fraser R. Flow (1000 m3/s)'].values
    sr = fdf['Skagit R. Flow (1000 m3/s)'].values
    this_yd = fdf.loc[TM, 'yearday']
    ax = fig.add_subplot(4,3,12)
    ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
    ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
    ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)
    # time markers
    ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(.9, .85, 'Columbia River', transform=ax.transAxes,
        color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .70, 'Fraser River', transform=ax.transAxes,
        color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .55, 'Skagit River', transform=ax.transAxes,
        color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-5,20)
    ax.set_axis_off()

    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = 'r'
    if tm.month in [4, 5, 6]:
        clist[1] = 'r'
    if tm.month in [7, 8, 9]:
        clist[2] = 'r'
    if tm.month in [10, 11, 12]:
        clist[3] = 'r'
    ax.text(0, 0, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=fs, horizontalalignment='left', style='italic')
    ax.text(.4, 0, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(.68, 0, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(1, 0, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=fs, horizontalalignment='right', style='italic')

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()

def P_superplot_chl(in_dict):
    # Plot phytoplankton maps and section, with forcing time-series.
    # Super clean design.  Updated to avoid need for tide data, which it
    # now just gets from the same mooring extraction it uses for wind.

    vn = 'phytoplankton'
    vlims = (0, 25) # full map
    vlims2 = (0, 25) # PS map
    vlims3 = (0, 25) # PS section
    cmap = 'Spectral_r'
    
    # get model fields
    ds = xr.open_dataset(in_dict['fn'])
    
    gtagex = str(in_dict['fn']).split('/')[-3]
    year_str = str(in_dict['fn']).split('/')[-2].split('.')[0][1:]

    # get forcing fields
    ffn = Ldir['LOo'] / 'extract' / gtagex / 'superplot' / ('forcing_' + gtagex + '_' + year_str + '.p')
    fdf = pd.read_pickle(ffn)
    fdf['yearday'] = fdf.index.dayofyear - 0.5 # .5 to 364.5

    # get section
    G, S, T = zrfun.get_basic_info(in_dict['fn'])
    # read in a section (or list of sections)
    tracks_path = Ldir['data'] / 'section_lines'
    tracks = ['Line_ps_main_v0.p']
    zdeep = -300
    xx = np.array([])
    yy = np.array([])
    for track in tracks:
        track_fn = tracks_path / track
        # get the track to interpolate onto
        pdict = pickle.load(open(track_fn, 'rb'))
        xx = np.concatenate((xx,pdict['lon_poly']))
        yy = np.concatenate((yy,pdict['lat_poly']))
    for ii in range(len(xx)-1):
        x0 = xx[ii]
        x1 = xx[ii+1]
        y0 = yy[ii]
        y1 = yy[ii+1]
        nn = 20
        if ii == 0:
            x = np.linspace(x0, x1, nn)
            y = np.linspace(y0,y1, nn)
        else:
            x = np.concatenate((x, np.linspace(x0, x1, nn)[1:]))
            y = np.concatenate((y, np.linspace(y0, y1, nn)[1:]))
    v2, v3, dist, idist0 = pfun.get_section(ds, vn, x, y, in_dict)

    # PLOTTING
    fig = plt.figure(figsize=(17,9))
    fs = 18 # fontsize

    # Full map
    ax = fig.add_subplot(131)
    lon = ds['lon_psi'].values
    lat = ds['lat_psi'].values
    v =ds[vn][0, -1, 1:-1, 1:-1].values
    fac=pinfo.fac_dict[vn]
    vv = fac * v
    vv[:, :6] = np.nan
    vv[:6, :] = np.nan
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims[0], vmax=vlims[1], cmap=cmap)
    pfun.add_coast(ax)
    ax.axis(pfun.get_aa(ds))
    pfun.dar(ax)
    ax.set_axis_off()
    # add a box for the subplot
    aa = [-123.5, -122.1, 47.03, 48.8]
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    # labels
    ax.text(.95, .07, 'LiveOcean\nPhytoplankton\n'+pinfo.units_dict[vn]+'\n'
        + datetime.strftime(T['dt'], '%Y'), fontsize=fs, color='k',
        transform=ax.transAxes, horizontalalignment='center',
        fontweight='bold')
    ax.text(.99,.97,'range\n'+ str(vlims), transform=ax.transAxes,
        va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    ax.text(.95, .03, datetime.strftime(T['dt'], '%Y.%m.%d'), fontsize=fs*.7, color='k',
        transform=ax.transAxes, horizontalalignment='center')

    # PS map
    ax = fig.add_subplot(132)
    cs = ax.pcolormesh(lon, lat, vv, vmin=vlims2[0], vmax=vlims2[1],
        cmap=cmap)
    #fig.colorbar(cs)
    pfun.add_coast(ax)
    ax.axis(aa)
    pfun.dar(ax)
    pfun.draw_box(ax, aa, color='c', alpha=.5, linewidth=5, inset=.01)
    ax.set_axis_off()
    # add section track
    sect_color = 'violet'
    n_ai = int(len(x)/6)
    n_tn = int(4.5*len(x)/7)
    ax.plot(x, y, linestyle='--', color='k', linewidth=2)
    ax.plot(x[n_ai], y[n_ai], marker='*', color=sect_color, markersize=14,
        markeredgecolor='k')
    ax.plot(x[n_tn], y[n_tn], marker='o', color=sect_color, markersize=10,
        markeredgecolor='k')
    # ax.text(.93,.97,'S range\n'+ str(vlims2), transform=ax.transAxes,
    #     va='top', ha='right', c='orange', size=.6*fs, weight='bold')
    

    # Section
    ax =  fig.add_subplot(433)
    ax.plot(dist, v2['zeta']+5, linestyle='--', color='k', linewidth=2)
    ax.plot(dist[n_ai], v2['zeta'][n_ai] + 5, marker='*', color=sect_color,
        markersize=14, markeredgecolor='k')
    ax.plot(dist[n_tn], v2['zeta'][n_tn] + 5, marker='o', color=sect_color,
        markersize=10, markeredgecolor='k')
    ax.set_xlim(dist.min(), dist.max())
    ax.set_ylim(zdeep, 25)
    sf = pinfo.fac_dict[vn] * v3['sectvarf']
    # plot section
    cs = ax.pcolormesh(v3['distf'], v3['zrf'], sf,
                       vmin=vlims3[0], vmax=vlims3[1], cmap=cmap)
    # ax.text(.99,.4,'S range\n'+ str(vlims3), transform=ax.transAxes,
    #     va='bottom', ha='right', c='orange', size=.6*fs, weight='bold')
                       
    #fig.colorbar(cs)
    # labels
    ax.text(0, 0, 'SECTION\nPuget Sound', fontsize=fs, color='b',
        transform=ax.transAxes)
    ax.set_axis_off()

    # get the day
    tm = T['dt'] # datetime
    TM = datetime(tm.year, tm.month, tm.day)
    # get yearday
    yearday = fdf['yearday'].values
    this_yd = fdf.loc[TM, 'yearday']

    # Tides
    alpha = .4
    ax = fig.add_subplot(436)
    ax.plot(yearday, fdf['RMS Tide Height (m)'].values, '-k',
        lw=3, alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM, 'RMS Tide Height (m)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(1, .05, 'NEAP TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    ax.text(1, .85, 'SPRING TIDES', transform=ax.transAxes,
        alpha=alpha, fontsize=fs, horizontalalignment='right')
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(0,1.5)
    ax.set_axis_off()

    # Wind
    alpha=.5
    ax = fig.add_subplot(439)
    w = fdf['8-day NS Wind Stress (Pa)'].values
    wp = w.copy()
    wp[w<0] = np.nan
    wm = w.copy()
    wm[w>0] = np.nan
    tt = np.arange(len(w))
    ax.fill_between(yearday, wp, y2=0*w, color='g', alpha=alpha)
    ax.fill_between(yearday, wm, y2=0*w, color='b', alpha=alpha)
    # time marker
    ax.plot(this_yd, fdf.loc[TM,'8-day NS Wind Stress (Pa)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(0, .85, 'DOWNWELLING WIND', transform=ax.transAxes,
        color='g', alpha=alpha, fontsize=fs)
    ax.text(0, .05, 'UPWELLING WIND', transform=ax.transAxes,
        color='b', alpha=alpha, fontsize=fs)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-.15, .25)
    ax.set_axis_off()

    # Rivers
    alpha = .6
    cr = fdf['Columbia R. Flow (1000 m3/s)'].values
    fr = fdf['Fraser R. Flow (1000 m3/s)'].values
    sr = fdf['Skagit R. Flow (1000 m3/s)'].values
    this_yd = fdf.loc[TM, 'yearday']
    ax = fig.add_subplot(4,3,12)
    ax.fill_between(yearday, cr, 0*yearday, color='orange', alpha=alpha)
    ax.fill_between(yearday, fr, 0*yearday, color='violet', alpha=alpha)
    ax.fill_between(yearday, sr, 0*yearday, color='brown', alpha=alpha)
    # time markers
    ax.plot(this_yd, fdf.loc[TM, 'Columbia R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Fraser R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    ax.plot(this_yd, fdf.loc[TM, 'Skagit R. Flow (1000 m3/s)'],
        marker='o', color='r', markersize=7)
    # labels
    ax.text(.9, .85, 'Columbia River', transform=ax.transAxes,
        color='orange', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .70, 'Fraser River', transform=ax.transAxes,
        color='violet', fontsize=fs, horizontalalignment='right', alpha=alpha)
    ax.text(.9, .55, 'Skagit River', transform=ax.transAxes,
        color='brown', fontsize=fs, horizontalalignment='right', alpha=alpha)
    # limits
    ax.set_xlim(0,365)
    ax.set_ylim(-5,20)
    ax.set_axis_off()

    # Time Axis
    clist = ['gray', 'gray', 'gray', 'gray']
    if tm.month in [1, 2, 3]:
        clist[0] = 'r'
    if tm.month in [4, 5, 6]:
        clist[1] = 'r'
    if tm.month in [7, 8, 9]:
        clist[2] = 'r'
    if tm.month in [10, 11, 12]:
        clist[3] = 'r'
    ax.text(0, 0, 'WINTER', transform=ax.transAxes, color=clist[0],
        fontsize=fs, horizontalalignment='left', style='italic')
    ax.text(.4, 0, 'SPRING', transform=ax.transAxes, color=clist[1],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(.68, 0, 'SUMMER', transform=ax.transAxes, color=clist[2],
        fontsize=fs, horizontalalignment='center', style='italic')
    ax.text(1, 0, 'FALL', transform=ax.transAxes, color=clist[3],
        fontsize=fs, horizontalalignment='right', style='italic')

    fig.tight_layout()

    # FINISH
    ds.close()
    if len(str(in_dict['fn_out'])) > 0:
        plt.savefig(in_dict['fn_out'])
        plt.close()
    else:
        plt.show()
