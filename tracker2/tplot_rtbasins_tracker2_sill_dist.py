"""
Plot results of a particle tracking experiment.
"""
from lo_tools import Lfun
from lo_tools import zfun
from lo_tools import plotting_functions as pfun
from lo_user_tools import llxyfun
Ldir = Lfun.Lstart()
import sys

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
from scipy.optimize import curve_fit
from scipy import stats
import tef_fun
import datetime

plt.close('all')
fig, axs = plt.subplots(1,5,figsize=(25,5))#,sharey=True)#,gridspec_kw={'height_ratios': [6,1]})
fig2, axs2 = plt.subplots(1,5,figsize=(25,5))#,sharey=True)#,gridspec_kw={'height_ratios': [6,1]})
fig3, ax3 = plt.subplots(1,1,figsize=(15,8))
# fig, ax = plt.subplots(1,1,figsize=(15,8))
inin_dist_plot=np.zeros(5)
outout_dist_plot=np.zeros(5)
inin_duration_plot=np.zeros(5)
outout_duration_plot=np.zeros(5)
silllens_plot=[5,10,20,40,80]

for i in range(5):
    # Choose an experiment and release to plot.
    # in_dir0 = Ldir['LOo'] / 'tracks'
    # exp_name = Lfun.choose_item(in_dir0, tag='', exclude_tag='.csv',
    #     itext='** Choose experiment from list **', last=False)
    # rel = Lfun.choose_item(in_dir0 / exp_name, tag='.nc', exclude_tag='grid',
    #     itext='** Choose item from list **', last=False)


    #for now just use 5,20,80 - can add 10 and 40 if it runs fast
    if i==0:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(45e3,0,45)
        linecolor = 'tab:red'
        silllenlabel = '5km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill5km_t0_xa0/sill5kmest_3d/release_2020.09.01.nc'
        tef_5km_fn = Ldir['LOo'] / 'extract/sill5km_t0_xa0/tef2/bulk_hourly_2020.09.01_2020.12.31' #this is for the Qprism timekeeper and shading
    elif i==1:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(50e3,0,45)
        linecolor = 'tab:orange'
        silllenlabel = '10km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill10km_t2_xa0/sill10kmest_3d/release_2020.09.01.nc'
    elif i==2:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(60e3,0,45)
        linecolor = 'tab:green'
        silllenlabel = '20km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill20kmdeep_t2_xa0/sill20kmdeepest_3d/release_2020.09.01.nc'
    elif i==3:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(80e3,0,45)
        linecolor = 'tab:blue'
        silllenlabel = '40km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill40km_t2_xa0/sill40kmest_3d/release_2020.09.01.nc'
    elif i==4:
        sillsea = llxyfun.x2lon(40e3,0,45)
        sillland = llxyfun.x2lon(120e3,0,45)
        linecolor = 'tab:purple'
        silllenlabel = '80km'
        fn = '/data1/ebroatch/LO_output/tracks2/sill80km_t2_xa0/sill80kmest_3d/release_2020.09.01.nc'
    
    print('\n'+silllenlabel+'\n')

    # get Datasets
    #fn = in_dir0 / exp_name / rel
    #fng = in_dir0 / exp_name / 'grid.nc'
    dsr = xr.open_dataset(fn, decode_times=False)
    #dsg = xr.open_dataset(fng)

    NT, NP = dsr.lon.shape

    # get a list of datetimes
    # ot_vec = dsr.ot.values
    # dt_list = [Lfun.modtime_to_datetime(ot) for ot in ot_vec]

    # # subsample output for plotting #SKIP SUBSAMPLING
    # npmax = 600 # max number of points to plot
    # step = max(1,int(np.floor(NP/npmax)))

    # lon = dsr.lon.values[:,::step]
    # lat = dsr.lat.values[:,::step]
    # lon = dsr.lon

    #longitude data
    lon_vals = dsr.lon.values #longitudes of the particles

    # z_vals = dsr.z.values #depths of the particles

    #time data
    time_hours = dsr.Time.values
    dsr.close()
    print('got lon_vals and time\n')

    # lon_start = lon_vals[np.newaxis, 0, :] #starting longitudes of the particles
    # lon_start_in = lon_start >= sillland #boolean array for particles starting in the inner basin
    # lon_start_insill = lon_start >= sillsea #boolean array for particles starting in the inner basin and sill
    # lon_start_out = (lon_start <= sillsea) & (lon_start >= 0) #boolean array for particles starting in the outer basin

    # z_start = z_vals[np.newaxis, 0, :] #starting depth of the particles
    # z_start_low = z_start < -50 #boolean array for particles starting below sill depth
    # z_start_up = z_start >= -50 #boolean array for particles starting above sill depth

    #find when particles are in different regions
    lon_in = lon_vals >= sillland #boolean array of particles in the inner basin over time
    lon_outocn = lon_vals <= sillsea #boolean array of particles in the outer basin or ocean over time
    lon_sill = (lon_vals > sillsea) & (lon_vals < sillland) #boolean array of particles on the sill over time

    #use the boolean array for particles on the sill to find on/off times
    sill_transition = np.diff(lon_sill.astype(int),axis=0) #must use astype int to differentiate between on (+1) and off (-1)
    #find the first and last on and off times in each column
    first_on = np.argmax(sill_transition,axis=0) #index of the first 1 in each column
    first_off = np.argmin(sill_transition,axis=0) #index of the first -1 in each column
    last_on = sill_transition.shape[0]-np.argmax(sill_transition[::-1,:],axis=0)-1 #flip array to get index of the last 1 in each column
    last_off = sill_transition.shape[0]-np.argmin(sill_transition[::-1,:],axis=0)-1 #flip array to get index of the last -1 in each column
    #we only want to use complete visits to the sill, so if the particle started or ended the run on the sill we want to remove that transition
    sill_transition_ends = np.zeros(sill_transition.shape)
    sill_transition_ends[:,:] = sill_transition #make a copy that we will alter the first and last on/off values if necessary
    #if the particle started on the sill, the first off (-1) will be before the first on (+1) and we will replace it with a 0
    sill_transition_ends[first_off,np.arange(sill_transition.shape[1])]=(first_on<first_off)*(-1) #keep the first off as -1 unless it comes before the first on
    #if the particle ended on the sill, the last on (1) will be after the last off (-1) and we will replace it with a 0
    sill_transition_ends[last_on,np.arange(sill_transition.shape[1])]=(last_on<last_off) #keep the last on as 1 unless it comes after the last off
    #because we use argmax and argmin, first and last indices will be set to zero or len-1 respectively if the on or off transition is not present
    #this causes problems if there is only one transition in the column, because the comparison of first_on<first_off and last_on<last_off will not detect the unpaired transition
    #therefore, if there is only one nonzero element in the column, set it to zero
    unpaired_par = np.where(np.count_nonzero(sill_transition_ends,axis=0)==1) #this gives the column/particle index of any unpaired single transition
    sill_transition_ends[:,unpaired_par]=0 #set those columns to all zero
    #now every column should have an equal number of on and off, and in each pair on comes before off
    #compared to the sill_times code, now we need to flatten everything so that we can use reduceat
    #this loses the particle index information
    #we will add an empty row to sill_transition_ends so that the indices of the labeled transitions correspond to the first on or off in the lon_vals
    sill_transition_ends_shift=np.concatenate((np.zeros((1,sill_transition_ends.shape[1])),sill_transition_ends),axis=0)
    sill_transition_ends_shift_flat=sill_transition_ends_shift.flatten(order='F')
    lon_vals_flat=lon_vals.flatten(order='F')
    #find the indices of on and off
    on_ind=np.where(sill_transition_ends_shift_flat==1)[0] #indices for on
    off_ind=np.where(sill_transition_ends_shift_flat==-1)[0] #indices for off
    #now we need to alternate the on and off indices to get a list of indices to make slices for reduceat
    #on should go first, and the list should be monotonically increasing
    reduce_inds=np.zeros(len(on_ind)+len(off_ind))
    reduce_inds[0::2]=on_ind
    reduce_inds[1::2]=off_ind
    #now we can find the maximum or minimum longitude in each slice interval using reduceat
    #every other interval will be the time off the sill so we will discard those
    max_lon_visit = np.maximum.reduceat(lon_vals_flat,reduce_inds.astype(int))
    min_lon_visit = np.minimum.reduceat(lon_vals_flat,reduce_inds.astype(int))
    max_lon_sill_visit = max_lon_visit[0::2]
    min_lon_sill_visit = min_lon_visit[0::2]
    #also get the durations
    durations=off_ind-on_ind
    #now we want to sort the visits by type using the region codes
    region_codes = (3*lon_in.astype(int))+(2*lon_sill.astype(int)) #NEW CODES: this gives 0 for outer basin and ocean, 2 for sill, and 3 for inner basin
    region_codes_flat = region_codes.flatten(order='F')
    on_from = region_codes_flat[on_ind-1] #need to subtract 1 because on_ind are the first time the particle is on the sill, we want the time before
    off_to = region_codes_flat[off_ind]
    #we are only interested in the max longitude for outer basin reflux and min for inner basin reflux
    inin_min_lon = np.where((on_from==3)&(off_to==3),min_lon_sill_visit,np.nan)
    outout_max_lon = np.where((on_from==0)&(off_to==0),max_lon_sill_visit,np.nan)
    inin_min_lon = inin_min_lon[~np.isnan(inin_min_lon)]
    outout_max_lon = outout_max_lon[~np.isnan(outout_max_lon)]
    print('average min lon inner basin reflux')
    print(np.mean(inin_min_lon))
    print('average max lon outer basin reflux:')
    print(np.mean(outout_max_lon))
    #convert the longitude into km from the end of the sill
    inin_dist_reached = (40+silllens_plot[i])-llxyfun.lon2x(inin_min_lon,0,45)/1e3
    outout_dist_reached = llxyfun.lon2x(outout_max_lon,0,45)/1e3 - 40
    print('average distance reached inner basin reflux')
    print(np.mean(inin_dist_reached))
    print('average distance reached outer basin reflux:')
    print(np.mean(outout_dist_reached))
    #also get the durations for the reflux visits
    inin_durations = np.where((on_from==3)&(off_to==3),durations,np.nan)
    outout_durations = np.where((on_from==0)&(off_to==0),durations,np.nan)
    inin_durations = inin_durations[~np.isnan(inin_durations)]
    outout_durations = outout_durations[~np.isnan(outout_durations)]
    #SAVE VALUES FOR PLOTTING LATER
    inin_dist_plot[i]=np.mean(inin_dist_reached)
    outout_dist_plot[i]=np.mean(outout_dist_reached)
    inin_duration_plot[i]=np.mean(inin_durations)
    outout_duration_plot[i]=np.mean(outout_durations)
    #PLOT HISTOGRAM

    #PLOT DISTANCE REACHED VS TIME ON SILL
    binmax=silllens_plot[i]
    binlist=np.arange(0,binmax+1,1)
    bincenters=np.arange(0.5,binmax+0.5,1)
    #now let's plot some histograms
    hist_inin,binedges_inin = np.histogram(inin_dist_reached,bins=binlist,density=True)
    hist_outout,binedges_outout = np.histogram(outout_dist_reached,bins=binlist,density=True)
    axs[i].plot(bincenters,hist_inin,lw=2,color='tab:pink',label='Inner basin reflux')
    axs[i].plot(bincenters,hist_outout,lw=2,color='tab:blue',label='Outer basin reflux')

    #make an additional plot all the histograms on one plot
    ax3.plot(bincenters,hist_inin,lw=2,color=linecolor,label=silllenlabel+' inner basin reflux')
    ax3.plot(bincenters,hist_outout,lw=2,color=linecolor,ls='--',label=silllenlabel+' outer basin reflux')

    #SCATTER PLOT OF DISTANCE VS DURATION
    axs2[i].scatter(inin_durations,inin_dist_reached,marker='o',color='tab:pink',alpha=0.3,label='Inner basin reflux')
    axs2[i].scatter(outout_durations,outout_dist_reached,marker='o',color='tab:blue',alpha=0.3,label='Outer basin reflux')
    durbinmax=np.max(durations)
    durbinlist=np.arange(-0.5,durbinmax+1.5,1)
    durbincenters=np.arange(0,durbinmax+1,1)
    inin_dist_durmean=stats.binned_statistic(inin_durations, inin_dist_reached, statistic='mean', bins=durbinlist, range=None)
    outout_dist_durmean=stats.binned_statistic(outout_durations, outout_dist_reached, statistic='mean', bins=durbinlist, range=None)
    axs2[i].plot(durbincenters,inin_dist_durmean.statistic,c=plt.cm.tab20(13),linewidth=2,ls='--',label='Inner basin reflux average')
    axs2[i].plot(durbincenters,outout_dist_durmean.statistic,c=plt.cm.tab20(19),linewidth=2,ls='--',label='Outer basin reflux average')

    # #now find the indices of on and off
    # on_ind=np.where(np.transpose(sill_transition_ends)==1)[::-1] #indices for on, sill_transition_ends[on_ind[0],on_ind[1]]=all 1's
    # off_ind=np.where(np.transpose(sill_transition_ends)==-1)[::-1] #indices for off, sill_transition_ends[on_ind[0],on_ind[1]]=all -1's
    # #the column (particle) index for both on and off will be identical, on_ind[1] is the same as off_ind[1]
    # #now find the times for on and off that will be compatible with the original lon_val array
    # on_times=on_ind[0]+1 #need to add 1 because on_ind was based on transitions from diff which is one element smaller, this gives the first time lon_vals in within sill
    # off_times=off_ind[0]+1
    # onoff_par=on_ind[1] #this gives the column/particle number that goes with the on and off times
    # #get the durations that particles spend on the sill
    # durations=off_times-on_times
    # print('average duration:')
    # print(np.mean(durations))
    # print('max duration:')
    # print(np.max(durations))
    # print('min duration:')
    # print(np.min(durations))
    # #debugging
    # # print('particle column mismatch:')
    # # print(np.count_nonzero(on_ind[1]!=off_ind[1]))
    # # print('number of first off to fix:')
    # # print(np.count_nonzero(first_off<first_on))
    # # print('number of first offs fixed:')
    # # print(np.count_nonzero(sill_transition==-1)-np.count_nonzero(sill_transition_ends==-1))
    # # print('number of last on to fix:')
    # # print(np.count_nonzero(last_off<last_on))
    # # print('number of last ons fixed:')
    # # print(np.count_nonzero(sill_transition==1)-np.count_nonzero(sill_transition_ends==1))
    
    # #split the durations by type of transit
    # #use a similar strategy to the returns fraction calculation
    # #make array with code number for each region
    # region_codes = (3*lon_in.astype(int))+(2*lon_sill.astype(int)) #NEW CODES: this gives 0 for outer basin and ocean, 2 for sill, and 3 for inner basin
    # print('got region codes\n')
    # #find transitions between regions
    # region_codes_transition = np.diff(region_codes,axis=0) #WITH NEW CODES this gives 0 staying in same region,-1 inner to sill,+1 sill to inner,+2 outer to sill,-2 sill to outer,+3 outer to inner direct,-3 inner to outer direct
    # #the region_codes array is the same size as lon_vals
    # #use the on_times and off_times to find which region the particles go to and from
    # on_from = region_codes[on_times-1,onoff_par] #need to subtract 1 because on_times are the first time the particle is on the sill, we want the time before
    # off_to = region_codes[off_times,onoff_par]
    # #now find the durations for only a single type of visit
    # inout_durations = np.where((on_from==3)&(off_to==0),durations,np.nan)
    # inin_durations = np.where((on_from==3)&(off_to==3),durations,np.nan)
    # outin_durations = np.where((on_from==0)&(off_to==3),durations,np.nan)
    # outout_durations = np.where((on_from==0)&(off_to==0),durations,np.nan)
    # #remove the nans
    # inout_durations = inout_durations[~np.isnan(inout_durations)]
    # inin_durations = inin_durations[~np.isnan(inin_durations)]
    # outin_durations = outin_durations[~np.isnan(outin_durations)]
    # outout_durations = outout_durations[~np.isnan(outout_durations)]
    # #these are the durations for each type of visit
    # #print some info
    # print('average duration in-out:')
    # print(np.mean(inout_durations))
    # print('average duration in-in:')
    # print(np.mean(inin_durations))
    # print('average duration out-in:')
    # print(np.mean(outin_durations))
    # print('average duration out-out:')
    # print(np.mean(outout_durations))

    # #up to this point we have ignored direct transits from inner-outer or outer-inner basins
    # #if the particle transits the sill in less than 1 hour, it may not be recorded on the sill in the hourly track file
    # #use the region codes to determine how many times this happens
    # #we can add these as zeros to the inout and outin duratoins
    # direct_inout_code = -3
    # direct_outin_code = +3
    # direct_inout_count = np.sum(region_codes_transition==direct_inout_code) #only matters for the 5km model
    # direct_outin_count = np.sum(region_codes_transition==direct_outin_code)
    # print('\nnumber of in->out direct: ')
    # print(direct_inout_count)
    # print('\nnumber of out->in direct: ')
    # print(direct_outin_count)
    # #append these <1h transits to the durations as zeros
    # inout_durations_full = np.concatenate((inout_durations,np.zeros(direct_inout_count)))
    # outin_durations_full = np.concatenate((outin_durations,np.zeros(direct_outin_count)))
    # print('average duration in-out including direct:')
    # print(np.mean(inout_durations_full))
    # print('average duration out-in including direct:')
    # print(np.mean(outin_durations_full))
    
    # #save some values for plotting later
    # inout_duration_plot[i]=np.mean(inout_durations_full)
    # inin_duration_plot[i]=np.mean(inin_durations)
    # outin_duration_plot[i]=np.mean(outin_durations_full)
    # outout_duration_plot[i]=np.mean(outout_durations)

    # binmax=np.max(durations)
    # binlist=np.arange(-0.5,binmax+1.5,1)
    # bincenters=np.arange(0,binmax+1,1)
    # #now let's plot some histograms
    # hist_inin,binedges_inin = np.histogram(inin_durations,bins=binlist,density=True)
    # hist_outout,binedges_outout = np.histogram(outout_durations,bins=binlist,density=True)
    # hist_inout,binedges_inout = np.histogram(inout_durations,bins=binlist,density=True)
    # hist_outin,binedges_outin = np.histogram(outin_durations,bins=binlist,density=True)
    # axs[i].plot(bincenters,hist_inin,lw=2,color='tab:pink',label='Inner basin reflux')
    # axs[i].plot(bincenters,hist_outout,lw=2,color='tab:blue',label='Outer basin reflux')
    # axs[i].plot(bincenters,hist_inout,lw=2,color=plt.cm.tab20(13),label='Efflux from inner basin')
    # axs[i].plot(bincenters,hist_outin,lw=2,color=plt.cm.tab20(19),label='Efflux from outer basin')
    # axs[i].hist(inin_durations,bins=binlist,histtype='step',color='tab:pink',label='Inner basin reflux')
    # axs[i].hist(outout_durations,bins=binlist,histtype='step',color='tab:blue',label='Outer basin reflux')
    # axs[i].hist(inout_durations,bins=binlist,histtype='step',color=plt.cm.tab20(13),label='Efflux from inner basin')
    # axs[i].hist(outin_durations,bins=binlist,histtype='step',color=plt.cm.tab20(19),label='Efflux from outer basin')

    # #make array with code number for each region
    # region_codes = lon_in.astype(int)+(2*lon_sill.astype(int)) #this gives 0 for outer basin and ocean, 2 for sill, and 1 for inner basin
    # print('got region codes\n')
    # #find transitions between regions
    # region_codes_transition = np.diff(region_codes,axis=0) #this gives 0 for staying in same region, +1 for inner to sill, -1 for sill to inner, +2 for outer to sill, -2 for sill to outer
    # #get the transitions all in a row with zeros removed
    # region_codes_transition_nan = np.where(region_codes_transition==0,np.nan,region_codes_transition) #change 0 to nan so that we can remove them and only look at consecutive transitions
    # a = (~np.isnan(region_codes_transition_nan)).argsort(0, kind='mergesort') #this gives the indices to sort the array with the nans first along the time axis, use mergesort to preserve order of other elements
    # region_codes_transition_consecutive = region_codes_transition_nan[a, np.arange(a.shape[1])[None,:]] #this should sort all the nans to the top of the column
    # print('got transition array\n')

    # #now we need to search for different patterns within the columns which indicate visits to the sill as efflux or reflux, and count them
    # #+1,-2 is in->sill->out (efflux)
    # #+1,-1 is in->sill->in (reflux to inner basin)
    # #+2,-1 is out->sill->in (efflux)
    # #+2,-2 is out->sill->out (reflux to outer basin)
    # pattern_inout = [1,-2]
    # pattern_inin = [1,-1]
    # pattern_outin = [2,-1]
    # pattern_outout = [2,-2]

    # inout_bool = (region_codes_transition_consecutive[:-1,:]==pattern_inout[0]) & (region_codes_transition_consecutive[1:,:]==pattern_inout[1]) #boolean arrays of where the patterns are found
    # inin_bool = (region_codes_transition_consecutive[:-1,:]==pattern_inin[0]) & (region_codes_transition_consecutive[1:,:]==pattern_inin[1])
    # outin_bool = (region_codes_transition_consecutive[:-1,:]==pattern_outin[0]) & (region_codes_transition_consecutive[1:,:]==pattern_outin[1])
    # outout_bool = (region_codes_transition_consecutive[:-1,:]==pattern_outout[0]) & (region_codes_transition_consecutive[1:,:]==pattern_outout[1])

    # inout_count = np.sum(inout_bool)
    # inin_count = np.sum(inin_bool)
    # outin_count = np.sum(outin_bool)
    # outout_count = np.sum(outout_bool)

    # print('got pattern counts\n')
    # print('\nin->sill->out: ')
    # print(inout_count)
    # print('\nin->sill->in: ')
    # print(inin_count)
    # print('\nout->sill->in: ')
    # print(outin_count)
    # print('\nout->sill->out: ')
    # print(outout_count)

    # #next, find the efflux reflux fractions
    # alpha_24 = inout_count/(inout_count+inin_count)
    # alpha_34 = inin_count/(inout_count+inin_count)
    # print('\ninner basin reflux (alpha_34): ')
    # print(alpha_34)
    # alpha_31 = outin_count/(outin_count+outout_count)
    # alpha_21 = outout_count/(outin_count+outout_count)
    # print('\nouter basin reflux (alpha_21): ')
    # print(alpha_21)   

    # alpha_24_plot[i]=alpha_24
    # alpha_34_plot[i]=alpha_34
    # alpha_31_plot[i]=alpha_31
    # alpha_21_plot[i]=alpha_21

    # lon_insill = lon_vals >= sillsea #boolean array of particles in the inner basin over time
    # lon_est = lon_vals >= 0 #boolean array of particles in the whole estuary over time
    # lon_out = (lon_vals <= sillsea) & (lon_vals >= 0) #boolean array of particles in the outer basin over time



    # z_low = z_vals < -50 #boolean array of particles below sill depth over time
    # z_up = z_vals >= -50 #boolean array of particles above sill depth over time

    #this includes returning particles
    # start_in_stay_in = lon_in * lon_start_in #particles in the inner basin that started in the inner basin
    # start_insill_stay_insill = lon_insill * lon_start_insill #particles in the inner basin and sill that started in the inner basin or sill #could change this to use lon_start_in?
    # start_est_stay_est = lon_est
    # start_out_stay_out = lon_out * lon_start_out #particles in the outer basin that started in the outer basin

    # #count number of returns to inner basin
    # ret_exit = np.zeros(lon_in.shape) 
    # ret_exit[1:,:] = np.diff(lon_in.astype(int), axis=0) #-1 for exit, +1 for return !!!important to cast lon_in as int or this won't work, it will count exits and returns as 1!!!
    # ret = np.where(ret_exit==1, 1, 0) #only keep the returns as +1 for the first hour that the particle is back inside the domain
    # ret_total = np.cumsum(ret, axis=0) #the number of returns a particle has made
    # ret_in_total = np.where(lon_in==0,np.nan,ret_total) #remove the particles that are currently outside
    # ret_in_total = np.where(lon_start_in==0,np.nan,ret_in_total) #remove the particles that started outside

    # ret_0 = ret_in_total==0
    # ret_1 = ret_in_total==1
    # ret_2 = ret_in_total==2
    # ret_3 = ret_in_total==3
    # ret_4plus = ret_in_total>=4

    # print('got data and boolean arrays\n')
    
    # ret_count_0 = np.nansum(ret_0,axis=1)
    # ret_count_1 = np.nansum(ret_1,axis=1)
    # ret_count_2 = np.nansum(ret_2,axis=1)
    # ret_count_3 = np.nansum(ret_3,axis=1)
    # ret_count_4plus = np.nansum(ret_4plus,axis=1)

    # ret_count_0_ta = zfun.lowpass(ret_count_0, f='godin')
    # ret_count_1_ta = zfun.lowpass(ret_count_1, f='godin')
    # ret_count_2_ta = zfun.lowpass(ret_count_2, f='godin')
    # ret_count_3_ta = zfun.lowpass(ret_count_3, f='godin')
    # ret_count_4plus_ta = zfun.lowpass(ret_count_4plus, f='godin')



    # #to count particles that have never left, find the time of first exit
    # tmax = time_hours[-1]
    # #inner basin
    # rt_strict_in = np.argmin(lon_in,axis=0).astype('float') #first time the particle is outside the inner basin
    # rt_strict_in = np.where(rt_strict_in==0, tmax+1, rt_strict_in) #replace 0 with tmax+1, this sets the particles that never leave or ones that were released outside to tmax+1
    # rt_strict_in = rt_strict_in * lon_start_in #this resets the particles that are not released in the inner basin to zero (necessary?)
    # rt_ind_in = rt_strict_in.astype(int)
    # exit_mask_in = np.zeros(lon_vals.shape, bool)
    # exit_mask_in[rt_ind_in[rt_ind_in<tmax],np.flatnonzero(rt_ind_in<tmax)] = True #this is a mask with the first time a particle exits set as True
    # strict_mask_in = np.logical_not(np.cumsum(exit_mask_in, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False
    # #inner basin + sill
    # rt_strict_insill = np.argmin(lon_insill,axis=0).astype('float')
    # rt_strict_insill = np.where(rt_strict_insill==0, tmax+1, rt_strict_insill) #replace 0 with tmax+1, this sets the particles that never leave or ones that were released outside to tmax+1
    # rt_strict_insill = rt_strict_insill * lon_start_insill #this resets the particles that are not released in the inner basin or sill to zero (necessary?)
    # rt_ind_insill = rt_strict_insill.astype(int)
    # exit_mask_insill = np.zeros(lon_vals.shape, bool)
    # exit_mask_insill[rt_ind_insill[rt_ind_insill<tmax],np.flatnonzero(rt_ind_insill<tmax)] = True #this is a mask with the first time a particle exits set as True
    # strict_mask_insill = np.logical_not(np.cumsum(exit_mask_insill, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False
    # #whole estuary
    # rt_strict_est = np.argmin(lon_est,axis=0).astype('float')
    # rt_strict_est = np.where(rt_strict_est==0, tmax+1, rt_strict_est) #replace 0 with tmax+1 (every particle is in the estuary at t=0, these are particles that never leave estuary)
    # rt_ind_est = rt_strict_est.astype(int)
    # exit_mask_est = np.zeros(lon_vals.shape, bool)
    # exit_mask_est[rt_ind_est[rt_ind_est<tmax],np.flatnonzero(rt_ind_est<tmax)] = True #this is a mask with the first time a particle exits set as True
    # strict_mask_est = np.logical_not(np.cumsum(exit_mask_est, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False


    # rt_strict = np.argmax(lon_vals<0,axis=0).astype('float')
    # rt_strict = np.where(rt_strict==0, tmax+1, rt_strict) #replace 0 with tmax+1 (every particle is in the estuary at t=0, these are particles that never leave estuary)
    # rt_ind = rt_strict.astype(int)
    # exit_mask = np.zeros(lon_vals.shape, bool)
    # exit_mask[np.flatnonzero(rt_ind<tmax), rt_ind[rt_ind<tmax]] = True #this is a mask with the first time a particle exits set as True
    # strict_mask = np.logical_not(np.cumsum(exit_mask, axis=0, dtype=bool)) #this is a mask with all times after a particle's first exit set as False

    # start_in_stay_in_strict = start_in_stay_in*strict_mask_in
    # start_insill_stay_insill_strict = start_insill_stay_insill*strict_mask_insill
    # start_est_stay_est_strict = start_est_stay_est*strict_mask_est
    # start_out_stay_out_strict = start_in_stay_in*strict_mask #OH NO NEED TO CHANGE THE LOGIC HERE #THIS WILL ONLY WORK FOR THE WHOLE ESTUARY SINCE THAT IS HOW THE RT IS DEFINED!!!

    # start_inlow_stay_inlow = lon_in * z_low * lon_start_in * z_start_low #particles that started in the inner basin below sill depth and stayed there
    # start_inup_stay_inup = lon_in * z_up * lon_start_in * z_start_up #particles that started in the inner basin above sill depth and stayed there
    # start_outlow_stay_outlow = lon_out * z_low * lon_start_out * z_start_low #particles that started in the outer basin below sill depth and stayed there
    # start_outup_stay_outup = lon_out * z_up * lon_start_out * z_start_up #particles that started in the outer basin above sill depth and stayed there

    # lon_in_lower = lon_in * lon_in[np.newaxis, 0, :] * z_start_lower #particles in the inner basin that started in the inner basin below sill depth
    # lon_in_upper = lon_in * lon_in[np.newaxis, 0, :] * z_start_upper #particles in the inner basin that started in the inner basin above sill depth

    # lon_out_lower = lon_out * lon_out[np.newaxis, 0, :] * z_start_lower #particles in the outer basin that started in the outer basin below sill depth
    # lon_out_upper = lon_out * lon_out[np.newaxis, 0, :] * z_start_upper #particles in the outer basin that started in the outer basin above sill depth


    # par_in_lower = np.sum(lon_in_lower,axis=1)
    # par_in_upper = np.sum(lon_in_upper,axis=1)
    # par_out_lower = np.sum(lon_out_lower,axis=1)
    # par_out_upper = np.sum(lon_out_upper,axis=1)

    # count_start_inlow_stay_in = np.sum(start_inlow_stay_in,axis=1) #total particles in the inner basin that started in the inner basin below sill depth
    # count_start_inup_stay_in = np.sum(start_inup_stay_in,axis=1) #total particles in the inner basin that started in the inner basin above sill depth
    # count_start_outlow_stay_out = np.sum(start_outlow_stay_out,axis=1) #total particles in the outer basin that started in the outer basin below sill depth
    # count_start_outup_stay_out = np.sum(start_outup_stay_out,axis=1) #total particles in the outer basin that started in the outer basin above sill depth

    # count_start_inlow_stay_inlow = np.sum(start_inlow_stay_inlow,axis=1) #total particles that started in the inner basin below sill depth and stayed there
    # count_start_inup_stay_inup = np.sum(start_inup_stay_inup,axis=1) #total particles that started in the inner basin above sill depth and stayed there
    # count_start_outlow_stay_outlow = np.sum(start_outlow_stay_outlow,axis=1) #total particles that started in the outer basin below sill depth and stayed there
    # count_start_outup_stay_outup = np.sum(start_outup_stay_outup,axis=1) #total particles that started in the outer basin above sill depth and stayed there

    # count_in = np.sum(start_in_stay_in,axis=1) #total particles in the inner basin that started in the inner basin
    # count_insill = np.sum(start_insill_stay_insill,axis=1) #total particles in the inner basin and sill that started in the inner basin or sill
    # count_est = np.sum(start_est_stay_est,axis=1) #total particles in the estuary

    # count_in_strict = np.sum(start_in_stay_in_strict,axis=1) #particles that have never left the inner basin
    # count_insill_strict = np.sum(start_insill_stay_insill_strict,axis=1) #particles that have never left the inner basin + sill
    # count_est_strict = np.sum(start_est_stay_est_strict,axis=1) #particles that have never left the estuary


    # print('got particle counts\n')
    # lon1 = dsr.lon.where((dsr.lon.sel(Time=0)<sillsea),drop=True) #THESE ARE THE OUTER BASIN PARTICLES
    # print('got lon1\n')
    # sys.stdout.flush()
    #SKIP SILL PARTICLES FOR NOW
    # lon2 = dsr.lon.where((dsr.lon.sel(Time=0)>=sillsea) & (dsr.lon.sel(Time=0)<sillland),drop=True)
    # print('got lon2\n')
    # # sys.stdout.flush()
    # lon3 = dsr.lon.where((dsr.lon.sel(Time=0)>=sillland),drop=True) #THESE ARE THE INNER BASIN PARTICLES
    # print('got lon3\n')
    # sys.stdout.flush()




    #FOR NOW, SKIP PARTICLES INITIALIZED ON SILL, AND ONLY PLOT PARTICLES REMAINING IN ORIGINAL BASIN
    # par1_ocn=(lon1<0).astype(int).sum(dim='Particle')
    # par1_out=((lon1>=0) & (lon1<sillsea)).astype(int).sum(dim='Particle') #THESE ARE THE OUTER BASIN PARTICLES STILL IN OUTER BASIN
    # par1_sill=((lon1>=sillsea) & (lon1<sillland)).astype(int).sum(dim='Particle')
    # par1_in=(lon1>=sillland).astype(int).sum(dim='Particle')
    # print('got par1\n')
    # sys.stdout.flush()

    # par2_ocn=(lon2<0).astype(int).sum(dim='Particle')
    # par2_out=((lon2>=0) & (lon2<sillsea)).astype(int).sum(dim='Particle')
    # par2_sill=((lon2>=sillsea) & (lon2<sillland)).astype(int).sum(dim='Particle')
    # par2_in=(lon2>=sillland).astype(int).sum(dim='Particle')
    # print('got par2\n')
    # sys.stdout.flush()

    # par3_ocn=(lon3<0).astype(int).sum(dim='Particle')
    # par3_out=((lon3>=0) & (lon3<sillsea)).astype(int).sum(dim='Particle')
    # par3_sill=((lon3>=sillsea) & (lon3<sillland)).astype(int).sum(dim='Particle')
    # par3_in=(lon3>=sillland).astype(int).sum(dim='Particle') #THESE ARE THE INNER BASIN PARTICLES STILL IN INNER BASIN
    # print('got par3\n')
    # sys.stdout.flush()

    # lonest = (lon>0)
    # lonest = lonest.astype(int)
    # partest = lonest.sum(dim='Particle')
    # partest_ta = zfun.lowpass(partest.values,f='godin',nanpad=True)[35:-35]
    # tplot = partest.Time.values[35:-35]


    # ax1.plot(par1_ocn.Time/24, zfun.lowpass((par1_ocn/par1_out.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:green', label='Ocean')
    # ax1.plot(par1_out.Time/24, zfun.lowpass((par1_out/par1_out.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:cyan', label='Outer basin')
    # ax1.plot(par1_out.Time/24, zfun.lowpass((par1_out/par1_out.sel(Time=0)).values*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY OUTER PARTICLES REMAINING
    # ax1.plot(par1_sill.Time/24, zfun.lowpass((par1_sill/par1_out.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:purple', label='Sill')
    # ax1.plot(par1_in.Time/24, zfun.lowpass((par1_in/par1_out.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:pink', label='Inner basin')

    #ax.set_ylim(20000,35000)


    # ax2.plot(par2_ocn.Time/24, zfun.lowpass((par2_ocn/par2_sill.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:green', label='Ocean')
    # ax2.plot(par2_out.Time/24, zfun.lowpass((par2_out/par2_sill.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:cyan', label='Outer basin')
    # ax2.plot(par2_sill.Time/24, zfun.lowpass((par2_sill/par2_sill.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:purple', label='Sill')
    # ax2.plot(par2_in.Time/24, zfun.lowpass((par2_in/par2_sill.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:pink', label='Inner basin')


    # ax3.plot(par3_ocn.Time/24, zfun.lowpass((par3_ocn/par3_in.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:green', label='Ocean')
    # ax3.plot(par3_out.Time/24, zfun.lowpass((par3_out/par3_in.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:cyan', label='Outer basin')
    # ax3.plot(par3_sill.Time/24, zfun.lowpass((par3_sill/par3_in.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:purple', label='Sill')
    # ax3.plot(par3_in.Time/24, zfun.lowpass((par3_in/par3_in.sel(Time=0)).values*100, f='godin'), linestyle=linst, color='tab:pink', label='Inner basin')
    # ax3.plot(par3_in.Time/24, zfun.lowpass((par3_in/par3_in.sel(Time=0)).values*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY INNER PARTICLES REMAINING
    

    # axs[0,0].plot(time_hours/24, zfun.lowpass((par_in_lower/par_in_lower[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY INNER PARTICLES REMAINING
    # axs[0,1].plot(time_hours/24, zfun.lowpass((par_in_upper/par_in_upper[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY INNER PARTICLES REMAINING
    # axs[1,0].plot(time_hours/24, zfun.lowpass((par_out_lower/par_out_lower[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY OUTER PARTICLES REMAINING
    # axs[1,1].plot(time_hours/24, zfun.lowpass((par_out_upper/par_out_upper[0])*100, f='godin'), color=linecolor, label=silllenlabel) #PLOT ONLY OUTER PARTICLES REMAINING

    # axs[0,0].plot(time_hours/24, zfun.lowpass(par_in_lower, f='godin'), color=linecolor, label=silllenlabel) #TRY WITH TOTAL PARTICLE COUNTS
    # axs[0,1].plot(time_hours/24, zfun.lowpass(par_in_upper, f='godin'), color=linecolor, label=silllenlabel) 
    # axs[1,0].plot(time_hours/24, zfun.lowpass(par_out_lower, f='godin'), color=linecolor, label=silllenlabel) 
    # axs[1,1].plot(time_hours/24, zfun.lowpass(par_out_upper, f='godin'), color=linecolor, label=silllenlabel) 

    # axs[0,0].plot(time_hours/24, zfun.lowpass(count_start_inlow_stay_in, f='godin'), color=linecolor, label=silllenlabel) #TRY WITH TOTAL PARTICLE COUNTS AND ADD STRICT LAYER SORTING
    # axs[0,0].plot(time_hours/24, zfun.lowpass(count_start_inlow_stay_inlow, f='godin'), '--', color=linecolor, label=silllenlabel)
    # axs[0,1].plot(time_hours/24, zfun.lowpass(count_start_inup_stay_in, f='godin'), color=linecolor, label=silllenlabel) 
    # axs[0,1].plot(time_hours/24, zfun.lowpass(count_start_inup_stay_inup, f='godin'), '--', color=linecolor, label=silllenlabel) 
    # axs[1,0].plot(time_hours/24, zfun.lowpass(count_start_outlow_stay_out, f='godin'), color=linecolor, label=silllenlabel) 
    # axs[1,0].plot(time_hours/24, zfun.lowpass(count_start_outlow_stay_outlow, f='godin'), '--', color=linecolor, label=silllenlabel) 
    # axs[1,1].plot(time_hours/24, zfun.lowpass(count_start_outup_stay_out, f='godin'), color=linecolor, label=silllenlabel+' in basin') 
    # axs[1,1].plot(time_hours/24, zfun.lowpass(count_start_outup_stay_outup, f='godin'), '--', color=linecolor, label=silllenlabel+' in layer')

    # #add shading #SKIP SHADING FOR NOW, COULD ADD BACK
    # if i==0:
    #     sect_name='b3'
    #     pad=36
    #     tef_df, vn_list, vec_list = tef_fun.get_two_layer(tef_5km_fn, sect_name)
    #     tef_df['Q_prism']=tef_df['qprism']/1000
    #     Qprism = tef_df['Q_prism'].loc['2020-09-04':'2020-12-28'] #cut off extra pad because qprism uses two godin filters
    #     ot=tef_df.loc['2020-09-04':'2020-12-28'].index
    #     ot_hours_delta = (((ot - datetime.datetime(2020,9,1,0,0,0)).total_seconds())/3600).to_numpy()
    #     axs[0,0].set_ylim(0,42000)
    #     axs[0,1].set_ylim(0,21000)
    #     axs[0,2].set_ylim(0,21000)
    #     axs[1,0].plot(ot_hours_delta/24,Qprism.to_numpy(), color='tab:gray', linewidth=2) #cut off the weird ends
    #     axs[1,1].plot(ot_hours_delta/24,Qprism.to_numpy(), color='tab:gray', linewidth=2)
    #     axs[1,2].plot(ot_hours_delta/24,Qprism.to_numpy(), color='tab:gray', linewidth=2)
    #     axs[1,0].set_ylabel('$Q_{prism}$ (5km)\n$[10^{3}\ m^{3}s^{-1}]$')
    #     axs[1,0].set_yticks(ticks=[20,50,80])
    #     axs[1,1].set_yticks(ticks=[20,50,80])
    #     axs[1,2].set_yticks(ticks=[20,50,80])
    #     axs[1,0].set_ylim(20,80)
    #     axs[1,1].set_ylim(20,80)
    #     axs[1,2].set_ylim(20,80)
    #     # ax0.set_xlim(pd.Timestamp('2020-09-01'), pd.Timestamp('2020-12-31'))
    #     snmid=(np.max(Qprism)+np.min(Qprism))/2
    #     snbg=np.where(Qprism.to_numpy()>snmid, 1, 0)
    #     axs[0,0].pcolor(ot_hours_delta/24, axs[0,0].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True) #slight change to the shading
    #     axs[0,1].pcolor(ot_hours_delta/24, axs[0,1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
    #     axs[0,2].pcolor(ot_hours_delta/24, axs[0,1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
    #     axs[1,0].pcolor(ot_hours_delta/24, axs[1,0].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
    #     axs[1,1].pcolor(ot_hours_delta/24, axs[1,1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
    #     axs[1,2].pcolor(ot_hours_delta/24, axs[0,1].get_ylim(), np.tile(snbg,(2,1)), cmap='Greys', vmin=-0.5, vmax=1.75, alpha=0.3, linewidth=0, antialiased=True)
    #     axs[1,0].grid(True)
    #     axs[1,1].grid(True)
    #     axs[1,2].grid(True)

    # axs[0,0].plot(time_hours/24, zfun.lowpass(count_est, f='godin'), color=linecolor, label=silllenlabel+' total')
    # axs[0,0].plot(time_hours/24, zfun.lowpass(count_est_strict, f='godin'), '--', color=linecolor, label=silllenlabel+' no return')
    # axs[0,1].plot(time_hours/24, zfun.lowpass(count_insill, f='godin'), color=linecolor, label=silllenlabel+' total') 
    # axs[0,1].plot(time_hours/24, zfun.lowpass(count_insill_strict, f='godin'), '--', color=linecolor, label=silllenlabel+' no return') 
    # axs[0,2].plot(time_hours/24, zfun.lowpass(count_in, f='godin'), color=linecolor, label=silllenlabel+' total') 
    # axs[0,2].plot(time_hours/24, zfun.lowpass(count_in_strict, f='godin'), '--', color=linecolor, label=silllenlabel+' no return') 

    # axs[i].stackplot(time_hours/24,ret_count_0_ta,ret_count_1_ta,ret_count_2_ta,ret_count_3_ta,ret_count_4plus_ta,colors=['tab:grey',plt.cm.tab20b(2),'tab:cyan','tab:olive','tab:pink'],labels=['0','1','2','3','4+'])

    #could try with total number of particles and/or double axis
    
    # axs[0,0].plot(time_hours/24, (par_out/par_out[0])*100, color=linecolor, label=silllenlabel) #TRY WITH NO FILTERING
    # axs[0,1].plot(time_hours/24, (par_in/par_in[0])*100, color=linecolor, label=silllenlabel)
    # print('plotted\n')
    sys.stdout.flush()
    
    #dsg.close()

# Add details to histogram plot
axs[0].set_title('5 km',color='tab:red')
axs[1].set_title('10 km',color='tab:orange')
axs[2].set_title('20 km',color='tab:green')
axs[3].set_title('40 km',color='tab:blue')
axs[4].set_title('80 km',color='tab:purple')
# axs[0].set_ylabel('# of particles')
axs[0].set_ylabel('fraction of particles')
axs[0].set_xlabel('Distance reached [km]')
axs[1].set_xlabel('Distance reached [km]')
axs[2].set_xlabel('Distance reached [km]')
axs[3].set_xlabel('Distance reached [km]')
axs[4].set_xlabel('Distance reached [km]')
axs[0].grid(True)
axs[1].grid(True)
axs[2].grid(True)
axs[3].grid(True)
axs[4].grid(True)
axs[0].set_xlim(0,5)
axs[1].set_xlim(0,10)
axs[2].set_xlim(0,20)
axs[3].set_xlim(0,40)
axs[4].set_xlim(0,80)
axs[0].set_ylim(bottom=0)
axs[1].set_ylim(bottom=0)
axs[2].set_ylim(bottom=0)
axs[3].set_ylim(bottom=0)
axs[4].set_ylim(bottom=0)
axs[4].legend()
fig.suptitle('Histograms of particle distances reached along sill')
fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_tracker2_sill_dist_hist.png' #UNCOMMENT TO PLOT
fig.savefig(fn_fig)

# Add details to single histogram plot
ax3.set_title('Histograms of particle distances reached along sill')
# axs[0].set_ylabel('# of particles')
ax3.set_ylabel('fraction of particles')
ax3.set_xlabel('Distance reached [km]')
ax3.grid(True)
ax3.set_xlim(0,80)
ax3.set_ylim(bottom=0)
ax3.legend()
fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_tracker2_sill_dist_hist_alt.png' #UNCOMMENT TO PLOT
fig3.savefig(fn_fig)

# Add details to scatter plot
axs2[0].set_title('5 km',color='tab:red')
axs2[1].set_title('10 km',color='tab:orange')
axs2[2].set_title('20 km',color='tab:green')
axs2[3].set_title('40 km',color='tab:blue')
axs2[4].set_title('80 km',color='tab:purple')
# axs[0].set_ylabel('# of particles')
axs2[0].set_ylabel('Distance reached [km]')
axs2[0].set_xlabel('Duration [h]')
axs2[1].set_xlabel('Duration [h]')
axs2[2].set_xlabel('Duration [h]')
axs2[3].set_xlabel('Duration [h]')
axs2[4].set_xlabel('Duration [h]')
axs2[0].grid(True)
axs2[1].grid(True)
axs2[2].grid(True)
axs2[3].grid(True)
axs2[4].grid(True)
axs2[0].set_xlim(0,14) #might need to adjust these
axs2[1].set_xlim(0,50)
axs2[2].set_xlim(0,100)
axs2[3].set_xlim(0,300)
axs2[4].set_xlim(0,400)
# axs2[0].set_ylim(0,5)
# axs2[1].set_ylim(0,10)
# axs2[2].set_ylim(0,20)
# axs2[3].set_ylim(0,40)
# axs2[4].set_ylim(0,80)
axs2[0].set_ylim(0,10) #just to check edges
axs2[1].set_ylim(0,15)
axs2[2].set_ylim(0,25)
axs2[3].set_ylim(0,45)
axs2[4].set_ylim(0,85)
axs2[4].legend()
fig2.suptitle('Particle distance reached by duration spent on sill')
fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_tracker2_sill_dist_scatter.png' #UNCOMMENT TO PLOT
fig2.savefig(fn_fig)

plt.close('all')

#ANOTHER FIGURE WITH THE AVERAGES ALL IN ONE PLOT
fig, ax = plt.subplots(1,1,figsize=(15,8))
axtwin = ax.twinx()
ax.plot(silllens_plot,inin_dist_plot,marker='o',c='tab:pink',ls='-',label=r'Distance reached (inner basin reflux)')
ax.plot(silllens_plot,outout_dist_plot,marker='o',c='tab:cyan',ls='-',label=r'Distance reached  (outer basin reflux)')
axtwin.plot(silllens_plot,inin_duration_plot,marker='^',c=plt.cm.tab20(13),ls='--',label=r'Duration (inner basin reflux)')
axtwin.plot(silllens_plot,outout_duration_plot,marker='^',c=plt.cm.tab20(19),ls='--',label=r'Duration (outer basin reflux)')
ax.set_xlabel('Sill length [km]')
ax.set_ylabel('Average distance reached along sill [km]')
axtwin.set_ylabel('Average duration spent on sill [h]')
ax.set_xlim(0,80)
ax.set_ylim(0,12)
axtwin.set_ylim(0,24)
ax.set_title('Extent along sill reached by refluxed particles')
ax.grid(True)
lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = axtwin.get_legend_handles_labels()
axtwin.legend(lines + lines2, labels + labels2, loc='lower right')
fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_tracker2_sill_dist_avg.png' #UNCOMMENT TO PLOT
fig.savefig(fn_fig)
plt.close()

# ax.plot(silllens_plot,alpha_34_plot,marker='o',c='tab:pink',ls='-',label=r'Inner basin reflux $\alpha_{34}$')
# ax.plot(silllens_plot,alpha_21_plot,marker='o',c='tab:cyan',ls='-',label=r'Outer basin reflux $\alpha_{21}$')
# ax.plot(silllens_plot,alpha_24_plot,marker='o',c=plt.cm.tab20(13),ls='--',label=r'Efflux from inner basin $\alpha_{24}$')
# ax.plot(silllens_plot,alpha_31_plot,marker='o',c=plt.cm.tab20(19),ls='--',label=r'Efflux from outer basin $\alpha_{31}$')
#plt.show()
#pfun.end_plot()

# #PLOTTING - HISTOGRAMS
# fig, axs = plt.subplots(5,1,sharex=True)
# for j in range(5):
#     hour=j*180
#     axs[j].set_title('t='+str(hour)+'h')
#     axs[j].hist(dsr['lon'].sel(Time=hour),bins=20,alpha=0.5)
#     #axs[j].set_ylim(0, 30)

#plt.show()
# ax.set_xlabel('Sill length [km]')
# ax.set_ylabel('Efflux/reflux coefficients')
# ax.set_xlim(0,80)
# ax.set_ylim(0,1)
# ax.set_title('Efflux/reflux fractions from particle trajectories')
# ax.grid(True)
# ax.legend()
# axs[0,0].set_xlabel('Days')
# axs[0,0].set_ylabel('% of particles remaining in inner basin')
# axs[0,0].set_ylabel('Particles remaining in inner basin') #TRY WITH TOTAL PARTICLE COUNT
# axs[0,0].set_ylabel('Particles remaining in inner basin') #TRY WITH TOTAL PARTICLE COUNT
# axs[0,0].set_title('Released in inner basin below sill height')
# axs[0,0].grid(True)
# axs[0,0].set_xlim(0,120)
# # axs[0,0].set_ylim(0,100)
# axs[0,0].set_ylim(0,count_start_inlow_stay_in[0])


# # axs[0,1].set_xlabel('Days')
# # axs[0,1].set_ylabel('% of particles')
# axs[0,1].set_title('Released in inner basin above sill height')
# axs[0,1].grid(True)
# axs[0,1].set_xlim(0,120)
# # axs[0,1].set_ylim(0,100)
# axs[0,1].set_ylim(0,count_start_inup_stay_in[0])

# axs[1,0].set_xlabel('Days')
# # axs[1,0].set_ylabel('% of particles remaining in outer basin')
# # axs[1,0].set_ylabel('Particles remaining in outer basin') #TRY WITH TOTAL PARTICLE COUNT
# axs[1,0].set_ylabel('Particles remaining') #TRY WITH TOTAL PARTICLE COUNT
# axs[1,0].set_title('Released in outer basin below sill height')
# axs[1,0].grid(True)
# axs[1,0].set_xlim(0,120)
# # axs[1,0].set_ylim(0,100)
# axs[1,0].set_ylim(0,count_start_outlow_stay_out[0])

# axs[1,1].set_xlabel('Days')
# # axs[1,1].set_ylabel('% of particles')
# axs[1,1].set_title('Released in outer basin above sill height')
# axs[1,1].grid(True)
# axs[1,1].set_xlim(0,120)
# # axs[1,1].set_ylim(0,100)
# axs[1,1].set_ylim(0,count_start_outup_stay_out[0])
# axs[1,1].legend(loc='upper right')

# ax2.set_xlabel('Days')
# #ax1.set_ylabel('Number of particles')
# ax2.set_title('Particles released on sill')
# #ax2.legend(loc='best')
# ax2.grid(True)
# ax2.set_xlim(0,120)
# ax2.set_ylim(0,100)

# ax3.set_xlabel('Days')
# #ax3.set_ylabel('Number of particles')
# ax3.set_title('Particles released in inner basin (unfiltered)')
# #ax3.legend(loc='best')
# ax3.grid(True)
# ax3.set_xlim(0,120)
# ax3.set_ylim(0,100)
# ax3.legend(loc='upper right')

# axs[0].grid(True)
# axs[1].grid(True)
# axs[2].grid(True)
# axs[3].grid(True)
# axs[4].grid(True)
# # axs[0,].grid(True)
# # axs[0,2].grid(True)

# axs[0].set_xlabel('Days')
# axs[1].set_xlabel('Days')
# axs[2].set_xlabel('Days')
# axs[3].set_xlabel('Days')
# axs[4].set_xlabel('Days')
# # axs[1,1].set_xlabel('Days')
# # axs[1,2].set_xlabel('Days')

# axs[0].set_xlim(0,120)
# axs[1].set_xlim(0,120)
# axs[2].set_xlim(0,120)
# axs[3].set_xlim(0,120)
# axs[4].set_xlim(0,120)
# # axs[0,1].set_xlim(0,120)
# # axs[0,2].set_xlim(0,120)
# # axs[1,0].set_xlim(0,120)
# # axs[1,1].set_xlim(0,120)
# # axs[1,2].set_xlim(0,120)

# axs[0].set_ylim(0,18000)
# # axs[0,1].set_ylim(0,21000)
# # axs[0,2].set_ylim(0,21000)

# # axs[0,0].set_title('Whole estuary')
# # axs[0,1].set_title('Inner basin + sill')
# # axs[0,2].set_title('Inner basin')

# axs[0].set_title('5 km',color='tab:red')
# axs[1].set_title('10 km',color='tab:orange')
# axs[2].set_title('20 km',color='tab:green')
# axs[3].set_title('40 km',color='tab:blue')
# axs[4].set_title('80 km',color='tab:purple')

# axs[0].set_ylabel('Particles remaining')
# axs[0].legend(title='Number of returns')
# axs[0,2].legend(loc='upper right')
#legend
# handles, labels = axs[0,0].get_legend_handles_labels()
# handles_reorder = np.concatenate((handles[::2],handles[1::2]),axis=0)
# labels_reorder = np.concatenate((labels[::2],labels[1::2]),axis=0)
# axs[0,0].legend(handles_reorder,labels_reorder,loc='upper right',ncol=2)


# fn_fig = Ldir['LOo'] / 'plots' / 'tplot_rtbasins_tracker2_sill_times.png' #UNCOMMENT TO PLOT
# plt.savefig(fn_fig)
# plt.close()

#plt.show()

