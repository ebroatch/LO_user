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
parser.add_argument('-g','--gridname', default='sill3', type=str)
parser.add_argument('-ctag','--collection_tag', default='test', type=str)
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

if args.gridname == 'sill3':
    sn_list = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
    x1km=np.array([0,10,20,30,38,40,42,44,46,48,50,58,68,78,86])
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

elif args.gridname == 'sill5km':
    sn_list = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
    x1km=np.array([0,10,20,30,38,40,41.25,42.5,43.75,45,47,55,65,75,83])
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
    sn_list = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
    x1km=np.array([0,10,20,30,38,40,42.5,45,47.5,50,52,60,70,80,88])
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
    sn_list = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
    x1km=np.array([0,10,20,30,38,40,45,50,55,60,62,70,80,90,98])
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
    sn_list = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
    x1km=np.array([0,10,20,30,38,40,50,60,70,80,82,90,100,110,118])
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
    sn_list = ['a1','a2','a3','a4','a5','b1','b2','b3','b4','b5','c1','c2','c3','c4','c5']
    x1km=np.array([0,10,20,30,38,40,60,80,100,120,122,130,140,150,158])
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