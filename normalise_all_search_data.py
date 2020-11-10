# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 18:03:28 2017

@author: psjohnsj
"""

#%% IMPORT LIBRARIES AND SET WORKING DIR
import os
import numpy as np
import pandas as pd
from cili.util import *
from cili.cleanup import *
from cili.extract import extract_events

datadir = 'C:\\Users\\psjohnsj\\Desktop\\pywkdir\\exp3.2_freeSearchCs\\data\\'
os.chdir(datadir)
ss = [sub for sub in os.walk(datadir).next()[1]]
sample_rate = 1000.0

#%%

def normalise_time(ranges, trial_info, values=None):
    interpx = np.linspace(0, 1, 100) # 100 datapoints in the new series
    sub_df = pd.DataFrame()
    all_trials = ranges[values].unstack(level=[0,1])
    for i in ranges.unstack().index:
        trial_df = pd.DataFrame({'time':range(0,100,1)})
        for val in values:
            d = all_trials.loc[val,i]
            x1 = np.arange(0, d.index.max()+1, dtype=np.float)
            rescale = x1/len(x1)
            interpy = np.interp(interpx, rescale, d.values)
            trial_df[val] = interpy   
            trial_df['trial'] = i
        for ti in trial_info:
            trial_df[ti] = ranges.loc[i,0][ti]
        sub_df = sub_df.append(trial_df)
    #sub_df.set_index(['trial','time'], inplace=True)
    return sub_df
    
#%%    

range_store = pd.HDFStore('ranges_r.h5')
new_norm_ranges = pd.HdFStore('new_norm_ranges.h5')
normalized_ranges = pd.HDFStore('normalized_ranges_r.h5')
trial_info = ['ntargets','blockNum','CORRECT','subject','TRIAL_INDEX','RT']
values     = ['c_pup_l_change','pup_l_change','x_l','y_l','pup_l','c_pup_l']

subs = list(set([val[:-2] for val in ss]))
all_subs_df = pd.DataFrame()    
for df in range_store.keys():
    rs = range_store[key]
    rs.RT = pd.to_numeric(rs.RT)
    rs.reset_index(inplace=True)
    new = rs.groupby('event', as_index=False, group_keys=False)\
            .apply(lambda g: g[g.onset <= g.RT.max()])
    new.set_index(['event','onset'], inplace=True)
    normed = normalise_time(new, trial_info=trial_info, values=values)
    normed.set_index(['trial','time'], inplace=True)
    normed['subject'] = key[1:-2]
    normed['session'] = key[-1]
    normed['subsesh'] = key
    
    new_norm_ranges.put(key, normed)


new_norm_ranges.close()

#%%

new_norm_ranges = pd.HDFStore('new_norm_ranges.h5')

gdf = pd.DataFrame()

for key in new_norm_ranges.keys():
    d = new_norm_ranges[key]
    gdf = gdf.append(d)
    
gdf.loc[gdf.session=='2', 'TRIAL_INDEX'] += 76
gdf.reset_index(inplace=True)
gdf.to_csv('new_norm_data.csv')
gdf = gdf.loc[gdf.CORRECT==1]
ag_xy = gdf.groupby(['subject','time'], as_index=False).mean()
ag_xy 