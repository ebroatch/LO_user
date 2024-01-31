"""
User-specific code for pgrid.

You would edit the information to reflect whatever grid you are working on.

Copied from LO example version updated January 18, 2024. Older grids by EB are in LO_user/pgrid/gfun_user_OLD.py for reference
"""
import numpy as np
import pandas as pd
from scipy import special
from lo_tools import zfun, Lfun
from lo_user_tools import llxyfun as lxf

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

# This is the name of the grid that you are working on.
gridname = 'sill20kmdeep'

# default s-coordinate info (could override below)
s_dict = {'THETA_S': 4, 'THETA_B': 2, 'TCLINE': 10, 'N': 30,
        'VTRANSFORM': 2, 'VSTRETCHING': 4}

def make_initial_info(gridname=gridname):
    # Add an elif section for your grid.

    if gridname == 'test0':
        # A large grid, used as a test.
        dch = gfun.default_choices()
        aa = [-130, -122, 42, 52]
        res = 1000 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north','south','west']
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['srtm15plus','cascadia','nw_pacific','psdem',
               'ttp_patch','grays_harbor','willapa_bay']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'hc0':
        dch = gfun.default_choices()
        aa = [-123.2, -122.537, 47.3, 47.9]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['psdem']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'ai0':
        dch = gfun.default_choices()
        aa = [-122.82, -122.36, 47.758, 48.18]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north', 'south', 'east', 'west']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['psdem']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'so1':
        # South Sound, new version, 2023.04.12
        dch = gfun.default_choices()
        dch['z_offset'] = -1.3 # NAVD88 is 1.3 m below MSL at Seattle
        dch['excluded_rivers'] = ['skokomish']
        aa = [-123.13, -122.76, 47, 47.42]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['east']
        dch['nudging_days'] = (0.1, 1.0)
        
        # by setting a small min_depth were are planning to use
        # wetting and drying in ROMS, but maintaining positive depth
        # for all water cells
        dch['min_depth'] = 0.2 # meters (positive down)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['psdem']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'wgh2':
        # Willapa Bay and Grays Harbor nest
        dch = gfun.default_choices()
        aa = [-124.4,-123.7,46.35,47.1]
        res = 200 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        
        dch['z_offset'] = -2
        # The docs for nw_pacific say the vertical datum is "sea level" and for Willapa
        # Bay and Grays Harbor it is MLLW so to match
        # this we would use z_offset = 0 or -1, but the intention here is to make the z=0
        # level be higher up, so that we catch more of the intertidal when using
        # WET_DRY. This should be matched by a similar intervention to zeta in ocnN.
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['nudging_days'] = (0.1, 1.0)
        
        # by setting a small min_depth were are planning to use
        # WET_DRY in ROMS, but maintaining positive depth
        # for all water cells
        dch['min_depth'] = 0.2 # meters (positive down)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['nw_pacific','grays_harbor','willapa_bay']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'cas2k':
        # cas6 domain but with 2 km resolution
        dch = gfun.default_choices()
        aa = [-130, -122, 42, 52]
        res = 2000 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['nudging_days'] = (3.0, 60.0)
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        
        # Initialize bathymetry
        dch['t_list'] = ['srtm15plus','cascadia','nw_pacific','psdem']
        z = gfu.combine_bathy_from_sources(lon, lat, dch)
                
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'cas7':
        # based completely on cas6 except we carve out Agate Pass and
        # Swinomish Channel by hand. This is an example of working from an
        # existing grid.
        dch = gfun.default_choices()
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['nudging_days'] = (3.0, 60.0)
        Ldir = Lfun.Lstart()
        fn = Ldir['parent'] / 'LO_output' / 'pgrid' / 'cas6' / 'grid.nc'
        dch['maskfile_to_copy'] = fn
        dch['remove_islands'] = False
        dch['trim_grid'] = False

        import xarray as xr
        ds = xr.open_dataset(fn)
        z = -ds.h.values
        lon = ds.lon_rho.values
        lat = ds.lat_rho.values
        
        # The plan is to only run:
        # start_grid
        # make_mask
        # edit_mask
        # (don't run carve_rivers - just copy the file from cas6)
        # smooth_grid
        # make_extras
        # grid_to_LO
            
    elif gridname == 'ae0':
        # analytical model estuary
        dch = gfun.default_choices()
        lon_list = [-2, 0, 1, 2]
        x_res_list = [2500, 500, 500, 2500]
        lat_list = [43, 44.9, 45.1, 47]
        y_res_list = [2500, 500, 500, 2500]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list,
                                            lat_list, y_res_list)
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        dch['analytical'] = True
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['use_z_offset'] = False
        # tidy up dch
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']
        # make bathymetry by hand
        z = np.zeros(lon.shape)
        x, y = zfun.ll2xy(lon, lat, 0, 45)
        zshelf = x * 1e-3
        zestuary = -20 + 20*x/1e5 + 20/(1e4)*np.abs(y)
        z = zshelf
        mask = zestuary < z
        z[mask] = zestuary[mask]
        
        # create a river file by hand
        Ldir = Lfun.Lstart()
        dch['ctag'] = 'ae0_v0'
        ri_dir = Ldir['LOo'] / 'pre' / 'river1' / dch['ctag']
        Lfun.make_dir(ri_dir)
        gri_fn = ri_dir / 'river_info.csv'
        with open(gri_fn, 'w') as rf:
            rf.write('rname,usgs,ec,nws,ratio,depth,flow_units,temp_units\n')
            rf.write('creek0,,,,1.0,5.0,m3/s,degC\n')
        # and make a track for the river
        track_dir = ri_dir / 'tracks'
        Lfun.make_dir(track_dir)
        track_fn = track_dir / 'creek0.p'
        track_df = pd.DataFrame()
        NTR = 100
        if True:
            track_df['lon'] = np.linspace(0,4,NTR) # OK to go past edge of domain
            track_df['lat'] = 45*np.ones(NTR)
        else: # Debugging with N/S river channel
            track_df['lon'] = 0.25*np.ones(NTR)
            track_df['lat'] = np.linspace(45,44,NTR) # South to North river
        track_df.to_pickle(track_fn)
        # *** NOTE: TRACKS MUST GO FROM OCEAN TO LAND ***

    elif gridname == 'sill20km':
        # Erin new idealized model similar to Admiralty Inlet

        #dch
        dch = gfun.default_choices()
        dch['analytical'] = True
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['use_z_offset'] = False
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']
        dch['z_land'] = 0 

        # bathymetry parameters
        # fixed bathymetry parameters
        L_basin = 40000 #40km long basins
        W_max = 8000 # 8km wide
        D_max = 200 # 200m max depth
        TL_sill = 2000 # steepness of sills (transition length in m, approx 10% slope)
        TL_side = 2000 # steepness of sides and end (transition length in m, approx 10% slope)
        # variable bathymetry parameters
        L_sill = 20000 # length of sill (length of flat section)
        HR_sill = 0.75 # height ratio of sills
        CR_sill = 0.5 # constriction ratio at sills
        # calculate additional constants
        L_estuary = (2*L_basin)+L_sill
        x_sill = L_estuary/2
        stretch_sill = TL_sill/4
        shift_sill = L_sill/2
        stretch_side = TL_side/4

        # make coordinate grid
        # resolution constants
        res_est=320 #resolution inside estuary
        res_ocn=2500 #resolution outside estuary
        # x and lon in estuary
        x_list_est_1 = np.flip(np.arange(x_sill, -res_est, -res_est)[1:])
        x_list_est_2 = np.arange(x_sill, L_estuary+res_est, res_est)
        x_list_est = np.concatenate((x_list_est_1, x_list_est_2), axis=None)
        lon_list_est = lxf.x2lon(x_list_est,0,45)
        # y and lat in estuary
        y_list_est_1 = np.flip(np.arange(0, -12000, -res_est)[1:])
        y_list_est_2 = np.arange(0, 12000, res_est)
        y_list_est = np.concatenate((y_list_est_1, y_list_est_2), axis=None)
        lat_list_est = lxf.y2lat(y_list_est,45)
        # ocean lon and south lat
        lon_list_1 = [-lon_list_est[0], 4]
        x_res_list_1 = [res_est, res_ocn]
        lat_list_1 = [-lat_list_est[0], -43]
        y_res_list_1 = [res_est, res_ocn]
        Lon_vec_1, Lat_vec_1 = gfu.stretched_grid(lon_list_1, x_res_list_1, lat_list_1, y_res_list_1)
        Lon_vec_1=-np.flip(Lon_vec_1[1:])
        Lat_vec_1=-np.flip(Lat_vec_1[1:])
        # land lon and north lat
        lon_list_2 = [lon_list_est[-1], 2.5] #make domain a bit bigger to accomodate future longer sill models
        x_res_list_2 = [res_est, res_ocn]
        lat_list_2 = [lat_list_est[-1], 47]
        y_res_list_2 = [res_est, res_ocn]
        Lon_vec_2, Lat_vec_2 = gfu.stretched_grid(lon_list_2, x_res_list_2, lat_list_2, y_res_list_2)
        Lon_vec_2=Lon_vec_2[1:]
        Lat_vec_2=Lat_vec_2[1:]
        # combine lon and lat vectors
        Lon_vec = np.concatenate((Lon_vec_1,lon_list_est,Lon_vec_2), axis=None)
        Lat_vec = np.concatenate((Lat_vec_1,lat_list_est,Lat_vec_2), axis=None)
        # make grids
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        x, y = zfun.ll2xy(lon, lat, 0, 45)

        # make bathymetry by hand
        #end profile
        D_end = ((D_max/2)*(special.erf((-(x-L_estuary)/stretch_side)-2)))+(D_max/2)
        #constrictions
        W_constrict=W_max-(((CR_sill*W_max/2)*(-special.erf(np.abs((x-x_sill)/stretch_sill)-(2+shift_sill/stretch_sill))))+(CR_sill*W_max/2))
        W_bottom=W_constrict-(2*TL_side)
        shift_bottom=W_bottom/2
        #basin shape
        z_basin=((D_end/2)*(special.erf(np.abs(y/stretch_side)-(2+shift_bottom/stretch_side))))-(D_end/2)
        #sill shape
        z_sill=-(((HR_sill*D_max/2)*(special.erf(np.abs((x-x_sill)/stretch_sill)-(2+shift_sill/stretch_sill))))+(D_max*(1-HR_sill/2)))
        #shelf
        z_shelf = x * 1e-3
        z_shelf = np.where(x>=0,0,z_shelf)
        #combine surfaces
        z_estuary=np.maximum(z_basin,z_sill)
        z_estuary=np.where(z_estuary>=-0.25,0,z_estuary)
        z=np.minimum(z_estuary,z_shelf)
        
        # create a river file by hand
        Ldir = Lfun.Lstart()
        dch['ctag'] = 'sill20km_base' #need to check what to call this
        ri_dir = Ldir['LOo'] / 'pre' / 'river1' / dch['ctag']
        Lfun.make_dir(ri_dir)
        gri_fn = ri_dir / 'river_info.csv'
        with open(gri_fn, 'w') as rf:
            rf.write('rname,usgs,ec,nws,ratio,depth,flow_units,temp_units\n')
            rf.write('creek0,,,,1.0,5.0,m3/s,degC\n')
        # and make a track for the river
        track_dir = ri_dir / 'tracks'
        Lfun.make_dir(track_dir)
        track_fn = track_dir / 'creek0.p'
        track_df = pd.DataFrame()
        NTR = 100
        if True:
            track_df['lon'] = np.linspace(0,2.5,NTR) # OK to go past edge of domain
            track_df['lat'] = 45*np.ones(NTR)
        else: # Debugging with N/S river channel
            track_df['lon'] = 0.25*np.ones(NTR)
            track_df['lat'] = np.linspace(45,44,NTR) # South to North river
        track_df.to_pickle(track_fn)
        # *** NOTE: TRACKS MUST GO FROM OCEAN TO LAND ***

    elif gridname == 'sill20kmdeep':
        # Erin new idealized model similar to Admiralty Inlet
        # change minimum depth and land height definitions to try to eliminate blowups

        #dch
        dch = gfun.default_choices()
        dch['analytical'] = True
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['use_z_offset'] = False
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']
        dch['z_land'] = 0 
        dch['min_depth'] = 8 #double minimum depth to try to reduce blowups (could be from shallow cells with large tidal amplitude)

        # bathymetry parameters
        # fixed bathymetry parameters
        L_basin = 40000 #40km long basins
        W_max = 8000 # 8km wide
        D_max = 200 # 200m max depth
        TL_sill = 2000 # steepness of sills (transition length in m, approx 10% slope)
        TL_side = 2000 # steepness of sides and end (transition length in m, approx 10% slope)
        # variable bathymetry parameters
        L_sill = 20000 # length of sill (length of flat section)
        HR_sill = 0.75 # height ratio of sills
        CR_sill = 0.5 # constriction ratio at sills
        # calculate additional constants
        L_estuary = (2*L_basin)+L_sill
        x_sill = L_estuary/2
        stretch_sill = TL_sill/4
        shift_sill = L_sill/2
        stretch_side = TL_side/4

        # make coordinate grid
        # resolution constants
        res_est=320 #resolution inside estuary
        res_ocn=2500 #resolution outside estuary
        # x and lon in estuary
        x_list_est_1 = np.flip(np.arange(x_sill, -res_est, -res_est)[1:])
        x_list_est_2 = np.arange(x_sill, L_estuary+res_est, res_est)
        x_list_est = np.concatenate((x_list_est_1, x_list_est_2), axis=None)
        lon_list_est = lxf.x2lon(x_list_est,0,45)
        # y and lat in estuary
        y_list_est_1 = np.flip(np.arange(0, -12000, -res_est)[1:])
        y_list_est_2 = np.arange(0, 12000, res_est)
        y_list_est = np.concatenate((y_list_est_1, y_list_est_2), axis=None)
        lat_list_est = lxf.y2lat(y_list_est,45)
        # ocean lon and south lat
        lon_list_1 = [-lon_list_est[0], 4]
        x_res_list_1 = [res_est, res_ocn]
        lat_list_1 = [-lat_list_est[0], -43]
        y_res_list_1 = [res_est, res_ocn]
        Lon_vec_1, Lat_vec_1 = gfu.stretched_grid(lon_list_1, x_res_list_1, lat_list_1, y_res_list_1)
        Lon_vec_1=-np.flip(Lon_vec_1[1:])
        Lat_vec_1=-np.flip(Lat_vec_1[1:])
        # land lon and north lat
        lon_list_2 = [lon_list_est[-1], 2.5] #make domain a bit bigger to accomodate future longer sill models
        x_res_list_2 = [res_est, res_ocn]
        lat_list_2 = [lat_list_est[-1], 47]
        y_res_list_2 = [res_est, res_ocn]
        Lon_vec_2, Lat_vec_2 = gfu.stretched_grid(lon_list_2, x_res_list_2, lat_list_2, y_res_list_2)
        Lon_vec_2=Lon_vec_2[1:]
        Lat_vec_2=Lat_vec_2[1:]
        # combine lon and lat vectors
        Lon_vec = np.concatenate((Lon_vec_1,lon_list_est,Lon_vec_2), axis=None)
        Lat_vec = np.concatenate((Lat_vec_1,lat_list_est,Lat_vec_2), axis=None)
        # make grids
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        x, y = zfun.ll2xy(lon, lat, 0, 45)

        # make bathymetry by hand
        #end profile
        D_end = ((D_max/2)*(special.erf((-(x-L_estuary)/stretch_side)-2)))+(D_max/2)
        #constrictions
        W_constrict=W_max-(((CR_sill*W_max/2)*(-special.erf(np.abs((x-x_sill)/stretch_sill)-(2+shift_sill/stretch_sill))))+(CR_sill*W_max/2))
        W_bottom=W_constrict-(2*TL_side)
        shift_bottom=W_bottom/2
        #basin shape
        z_basin=((D_end/2)*(special.erf(np.abs(y/stretch_side)-(2+shift_bottom/stretch_side))))-(D_end/2)
        #sill shape
        z_sill=-(((HR_sill*D_max/2)*(special.erf(np.abs((x-x_sill)/stretch_sill)-(2+shift_sill/stretch_sill))))+(D_max*(1-HR_sill/2)))
        #shelf
        z_shelf = x * 1e-3
        z_shelf = np.where(x>=0,2,z_shelf) #make land flat 2m elevation past the coast
        #combine surfaces
        z_estuary=np.maximum(z_basin,z_sill)
        z_estuary=np.where(z_estuary>=-0.25,2,z_estuary) #make very shallow cells part of land and set land elevation to 2m
        z=np.minimum(z_estuary,z_shelf)
        
        # create a river file by hand
        Ldir = Lfun.Lstart()
        dch['ctag'] = 'sill20km_base' #need to check what to call this
        ri_dir = Ldir['LOo'] / 'pre' / 'river1' / dch['ctag']
        Lfun.make_dir(ri_dir)
        gri_fn = ri_dir / 'river_info.csv'
        with open(gri_fn, 'w') as rf:
            rf.write('rname,usgs,ec,nws,ratio,depth,flow_units,temp_units\n')
            rf.write('creek0,,,,1.0,5.0,m3/s,degC\n')
        # and make a track for the river
        track_dir = ri_dir / 'tracks'
        Lfun.make_dir(track_dir)
        track_fn = track_dir / 'creek0.p'
        track_df = pd.DataFrame()
        NTR = 100
        if True:
            track_df['lon'] = np.linspace(0,2.5,NTR) # OK to go past edge of domain
            track_df['lat'] = 45*np.ones(NTR)
        else: # Debugging with N/S river channel
            track_df['lon'] = 0.25*np.ones(NTR)
            track_df['lat'] = np.linspace(45,44,NTR) # South to North river
        track_df.to_pickle(track_fn)
        # *** NOTE: TRACKS MUST GO FROM OCEAN TO LAND ***
        
    else:
        print('Error from make_initial_info: unsupported gridname')
        return
        
    if dch['trim_grid']:
        # check for odd size of grid and trim if needed
        NR, NC = lon.shape
        if np.mod(NR,2) != 0:
            print('- trimming row from grid')
            lon = lon[:-1,:]
            lat = lat[:-1,:]
            z = z[:-1,:]
        if np.mod(NC,2) != 0:
            print('- trimming column from grid')
            lon = lon[:,:-1]
            lat = lat[:,:-1]
            z = z[:,:-1]
            
    return lon, lat, z, dch
    

