"""
User-specific code for pgrid.

You would edit the information to reflect whatever grid you are working on.
"""
import numpy as np
import pandas as pd
from scipy import special
from lo_tools import zfun, Lfun

import sys
from pathlib import Path
pth = Path(__file__).absolute().parent.parent.parent / 'LO' / 'pgrid'
if str(pth) not in sys.path:
    sys.path.append(str(pth))
import gfun_utility as gfu
import gfun

# This is the name of the grid that you are working on.
gridname = 'ae0'

# default s-coordinate info (could override below)
s_dict = {'THETA_S': 4, 'THETA_B': 2, 'TCLINE': 10, 'N': 30,
        'VTRANSFORM': 2, 'VSTRETCHING': 4}

# Set the gridname and tag to use when creating the Ldir paths.
# They are used for accessing the river tracks, which may be developed for one
# grid but reused in others.
if gridname in ['ai0','hc0', 'sal0', 'so0']:
    # these cases reuse (all or some of) the LiveOcean cas6 model rivers
    base_gridname = 'cas6'
    base_tag = 'v3'
elif gridname in ['ae0']:
    # for analytical cases we create the river info and track in
    # make_initial_info() below, but we still assign gridname and tag so
    # that they get saved in the right places
    base_gridname = 'ae0'
    base_tag = 'v0'
elif gridname in ['sill0']:
    # for analytical cases we create the river info and track in
    # make_initial_info() below, but we still assign gridname and tag so
    # that they get saved in the right places
    base_gridname = 'sill0'
    base_tag = 'v0'
elif gridname in ['sill1']:
    # for analytical cases we create the river info and track in
    # make_initial_info() below, but we still assign gridname and tag so
    # that they get saved in the right places
    base_gridname = 'sill1'
    base_tag = 'v0'

def make_initial_info(gridname=gridname):
    # Add an elif section for your grid.

    if gridname == 'sal0':
        # A Salish Sea grid, used as an example.
        dch = gfun.default_choices()
        aa = [-124, -122, 47, 49]
        res = 600 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['nudging_edges'] = ['north', 'west']
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'hc0':
        dch = gfun.default_choices()
        aa = [-123.2, -122.537, 47.3, 47.9]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['t_list'] = [dch['t_dir'] / 'psdem' / 'PS_27m.nc']
        dch['nudging_edges'] = ['north']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'ai0':
        dch = gfun.default_choices()
        aa = [-122.82, -122.36, 47.758, 48.18]
        res = 100 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['t_list'] = [dch['t_dir'] / 'psdem_10m' / 'PS_30m.nc']
        dch['nudging_edges'] = ['north', 'south', 'east', 'west']
        dch['nudging_days'] = (0.1, 1.0)
        
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
    elif gridname == 'so0':
        # South Sound
        dch = gfun.default_choices()
        dch['z_offset'] = -1.3 # NAVD88 is 1.3 m below MSL at Seattle
        dch['excluded_rivers'] = ['skokomish']
        aa = [-123.13, -122.76, 47, 47.42]
        res = 50 # target resolution (m)
        Lon_vec, Lat_vec = gfu.simple_grid(aa, res)
        dch['t_list'] = [dch['t_dir'] / 'srtm15' / 'topo15.nc',
                    dch['t_dir'] / 'psdem_10m' / 'PS_30m.nc']
        dch['nudging_edges'] = ['east']
        dch['nudging_days'] = (0.1, 1.0)
        # Make the rho grid.
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)
        # initialize the bathymetry array
        z = np.nan * lon
        # add bathymetry automatically from files
        for t_fn in dch['t_list']:
            print('\nOPENING BATHY FILE: ' + t_fn.name)
            tlon_vec, tlat_vec, tz = gfu.load_bathy_nc(t_fn)
            tlon, tlat = np.meshgrid(tlon_vec, tlat_vec)
            z_part = zfun.interp2(lon, lat, tlon, tlat, tz)
            # put good values of z_part in z
            z[~np.isnan(z_part)] = z_part[~np.isnan(z_part)]
        if dch['use_z_offset']:
            z = z + dch['z_offset']
            
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
        
        # create a river file
        Ldir = Lfun.Lstart(gridname=base_gridname, tag=base_tag)
        ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag']
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
        track_df['lon'] = np.linspace(0,4,NTR) # OK to go past edge of domain
        track_df['lat'] = 45*np.ones(NTR)
        track_df.to_pickle(track_fn)
        # NOTE: tracks go from ocean to land

    elif gridname == 'sill0':
        # analytical model estuary with sills
        # currently using parameters for the "base case"
        dch = gfun.default_choices()
        lon_list = [-4, 0, 2, 3, 4]
        x_res_list = [2500, 500, 500, 2500, 2500]
        lat_list = [43, 44.9, 45.1, 47]
        y_res_list = [2500, 500, 500, 2500]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list, lat_list, y_res_list)
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)

        dch['analytical'] = True
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['use_z_offset'] = False
        # tidy up dch
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']

        # fixed bathymetry parameters
        L_estuary = 160000 # 160km long
        W_estuary = 8000 # 8km wide
        D_estuary = 200 # 200m max depth
        sill_steep = 2000 # steepness of sills (transition length in m, approx 10% slope)
        side_steep = 2000 # steepness of sides and end (transition length in m, approx 10% slope)

        # variable bathymetry parameters
        sill_num = 3 # number of sills
        sill_ratio = 0.5 # height ratio of sills
        sill_length = 4000 # length of sill (length of flat section)
        constrict_ratio = 0 # constriction ratio at sills

        # calculate additional constants
        L_basin = L_estuary/(sill_num+1)
        sill_pos = np.linspace(L_basin,L_estuary-L_basin,3)
        sill_stretch = sill_steep/4
        sill_shift = sill_length/2
        side_stretch = side_steep/4

        # make bathymetry by hand
        #grid
        x, y = zfun.ll2xy(lon, lat, 0, 45)
        #end profile
        depth_end = ((D_estuary/2)*(special.erf((-(x-L_estuary)/side_stretch)-2)))+(D_estuary/2)
        #constrictions
        width_constrict=W_estuary*np.ones(Lon_vec.size)
        for i in range(sill_num):
            constrict=W_estuary-(((constrict_ratio*W_estuary/2)*(-special.erf(np.abs((x-sill_pos[i])/sill_stretch)-(2+sill_shift/sill_stretch))))+(constrict_ratio*W_estuary/2))
            width_constrict=np.minimum(width_constrict,constrict)
        bottom_width=width_constrict-(2*side_steep)
        bottom_shift=bottom_width/2
        #basin shape
        z_basin=((depth_end/2)*(special.erf(np.abs(y/side_stretch)-(2+bottom_shift/side_stretch))))-(depth_end/2)
        #sill shape
        z_sills=-D_estuary*np.ones(x.shape)
        for i in range(sill_num):
            sill=-(((sill_ratio*D_estuary/2)*(special.erf(np.abs((x-sill_pos[i])/sill_stretch)-(2+sill_shift/sill_stretch))))+(D_estuary*(1-sill_ratio/2)))
            z_sills=np.maximum(z_sills,sill)
        #shelf
        z_shelf = x * 1e-3
        #combine surfaces
        z_estuary=np.maximum(z_basin,z_sills)
        z=np.minimum(z_estuary,z_shelf)

        # create a river file
        Ldir = Lfun.Lstart(gridname=base_gridname, tag=base_tag)
        ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag']
        Lfun.make_dir(ri_dir)
        gri_fn = ri_dir / 'river_info.csv'
        with open(gri_fn, 'w') as rf:
            rf.write('rname,usgs,ec,nws,ratio,depth,flow_units,temp_units\n')
            rf.write('creek_sill0,,,,1.0,5.0,m3/s,degC\n')
        # and make a track for the river
        track_dir = ri_dir / 'tracks'
        Lfun.make_dir(track_dir)
        track_fn = track_dir / 'creek_sill0.p'
        track_df = pd.DataFrame()
        NTR = 100
        track_df['lon'] = np.linspace(0,4,NTR) # OK to go past edge of domain
        track_df['lat'] = 45*np.ones(NTR)
        track_df.to_pickle(track_fn)
        # NOTE: tracks go from ocean to land

    elif gridname == 'sill1':
        # analytical model estuary with sills
        # currently using parameters for the "base case"
        dch = gfun.default_choices()
        lon_list = [-4, 0, 2, 3, 4]
        x_res_list = [2500, 500, 500, 2500, 2500]
        lat_list = [43, 44.9, 45.1, 47]
        y_res_list = [2500, 500, 500, 2500]
        Lon_vec, Lat_vec = gfu.stretched_grid(lon_list, x_res_list, lat_list, y_res_list)
        lon, lat = np.meshgrid(Lon_vec, Lat_vec)

        dch['analytical'] = True
        dch['nudging_edges'] = ['north', 'south', 'west']
        dch['use_z_offset'] = False
        # tidy up dch
        dch['z_offset'] = 0.0
        dch['t_dir'] = 'BLANK'
        dch['t_list'] = ['BLANK']

        # fixed bathymetry parameters
        L_basin = 40000 #40km long basins
        W_max = 8000 # 8km wide
        D_max = 200 # 200m max depth
        TL_sill = 2000 # steepness of sills (transition length in m, approx 10% slope)
        TL_side = 2000 # steepness of sides and end (transition length in m, approx 10% slope)

        # variable bathymetry parameters
        L_sill = 4000 # length of sill (length of flat section)
        HR_sill = 0.75 # height ratio of sills
        CR_sill = 0 # constriction ratio at sills

        # calculate additional constants
        L_estuary = (2*L_basin)+L_sill
        x_sill = L_estuary/2
        stretch_sill = TL_sill/4
        shift_sill = L_sill/2
        stretch_side = TL_side/4

        # make bathymetry by hand
        #grid
        x, y = zfun.ll2xy(lon, lat, 0, 45)
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
        #combine surfaces
        z_estuary=np.maximum(z_basin,z_sill)
        z=np.minimum(z_estuary,z_shelf)

        # create a river file
        Ldir = Lfun.Lstart(gridname=base_gridname, tag=base_tag)
        ri_dir = Ldir['LOo'] / 'pre' / 'river' / Ldir['gtag']
        Lfun.make_dir(ri_dir)
        gri_fn = ri_dir / 'river_info.csv'
        with open(gri_fn, 'w') as rf:
            rf.write('rname,usgs,ec,nws,ratio,depth,flow_units,temp_units\n')
            rf.write('creek_sill1,,,,1.0,5.0,m3/s,degC\n')
        # and make a track for the river
        track_dir = ri_dir / 'tracks'
        Lfun.make_dir(track_dir)
        track_fn = track_dir / 'creek_sill1.p'
        track_df = pd.DataFrame()
        NTR = 100
        track_df['lon'] = np.linspace(0,4,NTR) # OK to go past edge of domain
        track_df['lat'] = 45*np.ones(NTR)
        track_df.to_pickle(track_fn)
        # NOTE: tracks go from ocean to land

    else:
        print('Error from make_initial_info: unsupported gridname')
        return
        
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
    

