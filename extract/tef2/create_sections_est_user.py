# THIS IS A VERSION OF create_sections_user WITH ONLY ONE SECTION (same as a1) AT THE MOUTH OF THE ESTUARY
# This will allow us to use the segment tools to get the mean salinity of the whole estuary
# For use in the salt and freshwater flushing times calculations

from lo_tools import Lfun
from lo_user_tools import llxyfun
import numpy as np
import pandas as pd
import xarray as xr
import sys

# get command line arguments
import argparse
parser = argparse.ArgumentParser()
# gridname and tag will create a collection folder
parser.add_argument('-g','--gridname', default='sill5km', type=str)
parser.add_argument('-ctag','--collection_tag', default='cest', type=str)
# set clobber to True to start a new collection even if it exists
parser.add_argument('-clobber', default=False, type=Lfun.boolean_string)
args = parser.parse_args()

# paths to save sections
Ldir = Lfun.Lstart(gridname=args.gridname)
out_name = 'sections_' + args.gridname + '_' + args.collection_tag
out_dir = Ldir['LOo'] / 'extract' / 'tef2' / out_name
if args.clobber:
    inp = input('Do you really want to clobber? y/n: ')
    if inp == 'y':
        Lfun.make_dir(out_dir, clean=True)
    else:
        sys.exit()
else:
    # make out_dir in case it does not exist yet
    Lfun.make_dir(out_dir)

if args.gridname == 'sill5km':
    sn_list = ['a1est']
    x1km=np.array([0])
    x1m=x1km*1e3
    x1=llxyfun.x2lon(x1m,0,45)
    y1=45.05
    y2=44.95
    for i in range(len(sn_list)):
        sn=sn_list[i]
        df = pd.DataFrame(columns=['x','y'])
        df.loc[:,'x'] = np.array([x1[i], x1[i]])
        df.loc[:,'y'] = np.array([y1, y2])
        df.to_pickle(out_dir / (sn + '.p'))

elif args.gridname == 'sill10km':
    sn_list = ['a1est']
    x1km=np.array([0])
    x1m=x1km*1e3
    x1=llxyfun.x2lon(x1m,0,45)
    y1=45.05
    y2=44.95
    for i in range(len(sn_list)):
        sn=sn_list[i]
        df = pd.DataFrame(columns=['x','y'])
        df.loc[:,'x'] = np.array([x1[i], x1[i]])
        df.loc[:,'y'] = np.array([y1, y2])
        df.to_pickle(out_dir / (sn + '.p'))

elif args.gridname == 'sill20kmdeep':
    sn_list = ['a1est']
    x1km=np.array([0])
    x1m=x1km*1e3
    x1=llxyfun.x2lon(x1m,0,45)
    y1=45.05
    y2=44.95
    for i in range(len(sn_list)):
        sn=sn_list[i]
        df = pd.DataFrame(columns=['x','y'])
        df.loc[:,'x'] = np.array([x1[i], x1[i]])
        df.loc[:,'y'] = np.array([y1, y2])
        df.to_pickle(out_dir / (sn + '.p'))

elif args.gridname == 'sill40km':
    sn_list = ['a1est']
    x1km=np.array([0])
    x1m=x1km*1e3
    x1=llxyfun.x2lon(x1m,0,45)
    y1=45.05
    y2=44.95
    for i in range(len(sn_list)):
        sn=sn_list[i]
        df = pd.DataFrame(columns=['x','y'])
        df.loc[:,'x'] = np.array([x1[i], x1[i]])
        df.loc[:,'y'] = np.array([y1, y2])
        df.to_pickle(out_dir / (sn + '.p'))

elif args.gridname == 'sill80km':
    sn_list = ['a1est']
    x1km=np.array([0])
    x1m=x1km*1e3
    x1=llxyfun.x2lon(x1m,0,45)
    y1=45.05
    y2=44.95
    for i in range(len(sn_list)):
        sn=sn_list[i]
        df = pd.DataFrame(columns=['x','y'])
        df.loc[:,'x'] = np.array([x1[i], x1[i]])
        df.loc[:,'y'] = np.array([y1, y2])
        df.to_pickle(out_dir / (sn + '.p'))

else:
    print('Error: unsupported gridname')