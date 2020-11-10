# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 13:00:09 2016

@author: Joel
"""

import numpy as np
import pandas as pd
#%% FUNCTIONS

def events_from_message_report(msg_rep, sub, message):
    
    sub_msgs = msg_rep[(msg_rep.subsesh==sub) &\
                       (msg_rep.CURRENT_MSG_TEXT==message)]
    idx = sub_msgs['TRIAL_START_TIME'] + sub_msgs['CURRENT_MSG_TIME']
    sub_msgs.index = idx
    
    return sub_msgs

def print_info(events, samples, ranges):
    
    print '....................'
    print 'EBLINKs            : ' + str(len(events.EBLINK))
    print 'EFIXs              : ' + str(len(events.EFIX))
    print '%_interp_overall   : ' + str(samples.interpolated.value_counts(normalize=True)[1]*100) 
    print '% interp_in_trials : ' + str(ranges.interpolated.value_counts(normalize=True)[1]*100)


    
    
def normalise_time(ranges, trial_info):

    all_trials = ranges['c_pup_l_change'].unstack()
    pup_dict={}
    for i in ranges.unstack().index:
        offset = int(ranges.loc[i, 'RT'][0])
        pup = all_trials.loc[i, 0:offset-1]
        x1 = np.arange(0, offset, dtype=np.float)
        rescale = x1/len(x1)
        interpx = np.linspace(0, 1, 100) # 100 datapoints in the new series
        interpy = np.interp(interpx, rescale, pup)
        pup_dict[i] = interpy
    info = ranges[trial_info]
    info.reset_index(inplace=True)
    info = info[info.onset < 100]
    info.set_index(['event','onset'], inplace=True)
    norm_onset = [range(100) for val in range(0,76) ]
    norm_onset = [val for sub in norm_onset for val in sub] 
    pup_df = pd.DataFrame(pup_dict)
    pup_df = pd.melt(pup_df, var_name='event', value_name='c_pup_l_change')
    pup_df['onset'] = norm_onset
    pup_df.set_index(['event','onset'], inplace=True)
    pup_df = pup_df.join(info)
    
    return pup_df

def get_fix_aligned(samples, events, trial_info):
    
    ranges        = extract_events(samples, 
                                   events,
                                   offset=-500,
                                   duration=2500,
                                   units='samples',
                                   borrow_attributes=trial_info)
    baselines     = extract_events(samples, 
                                   events,
                                   offset=-50,
                                   duration=100,
                                   units='samples',
                                   borrow_attributes=trial_info).mean(level=0)  
    
    ranges['pup_l_change']   = (((ranges.pup_l - baselines.pup_l) / baselines.pup_l) * 100).values # % signal change from baseline
    ranges['c_pup_l_change'] = (((ranges.c_pup_l - baselines.c_pup_l) / baselines.c_pup_l) * 100).values # % signal change from baseline
    ranges.reset_index(inplace=True)    
    
    return ranges
    
