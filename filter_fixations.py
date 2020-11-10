# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 10:55:50 2016

@author: Joel
"""
#%%
import os
import pandas as pd
import numpy as np
datadir = 'C:\\Users\\Joel\\Desktop\\pywkdir\\vissearchCs\\data\\'
os.chdir(datadir)

#%% FILTER FIXATIONS
#fixs = pd.read_table('fixation_report3.txt', low_memory=False)
#fixs.rename(columns={'RECORDING_SESSION_LABEL':'subsesh'}, inplace=True)
#fixs['subject'] = fixs.subsesh.str[:-2]
#fixs['session'] = fixs.subsesh.str[-1]
#fixs['timestamp'] = fixs.CURRENT_FIX_START + fixs.TRIAL_START_TIME
#fixs.to_csv('fixation_report3.csv')

fixs = pd.read_table('fixation_report3.csv', sep=',', na_values='.')
fixs['refix'] = 0
fixs.loc[fixs.CURRENT_FIX_REFIX_INTEREST_AREA != 0, 'refix'] = 1
fixs['fix_idx'] = fixs.index
fixs.index = fixs.timestamp.values
fixs = fixs[fixs.blockNum != 1] # remove practice trials
fixs = fixs[fixs.TRIAL_INDEX <= 76]
fixs = fixs[fixs.ntargets.notnull()]
target_msgs = ['TARGET_1','TARGET_2','TARGET_3']

# find out how many of the targets were fixated
fixs['ntargets_fixed'] = 0
for key, df in fixs.groupby(['subsesh','TRIAL_INDEX']):
    ntar = int(df.ntargets.unique()[0])
    ntf = 0
    for i in range(ntar):
        if target_msgs[i] in df.CURRENT_FIX_INTEREST_AREA_LABEL.values:
            ntf += 1
        fixs.loc[(fixs.subsesh==key[0])&(fixs.TRIAL_INDEX==key[1]), 'ntargets_fixed'] = ntf    
    
# get target fixations
tarfixs = fixs[(fixs.CURRENT_FIX_NEAREST_INTEREST_AREA_LABEL.isin(target_msgs)) & \
                (fixs.CURRENT_FIX_NEAREST_INTEREST_AREA_DISTANCE <= 1.50)]

# use only 1st fixations in interest area and those which happen within a certain timeframe               
tarfixs = tarfixs[tarfixs.CURRENT_FIX_RUN_INDEX == 1]
tarfixs = tarfixs[tarfixs.CURRENT_FIX_DURATION >= 120]
#tarfixs = tarfixs[tarfixs.CURRENT_FIX_START > 1000] # can do this later
#tarfixs = tarfixs[tarfixs.CURRENT_FIX_START  < tarfixs.RT - 1500] # can do this later

# trim any above mean * 3sd
tar_sd_dur = tarfixs.CURRENT_FIX_DURATION.std()
tar_mean_dur = tarfixs.CURRENT_FIX_DURATION.mean()
tarfixs = tarfixs[tarfixs.CURRENT_FIX_DURATION <= tar_mean_dur + tar_sd_dur*3]
tar_max_dur = tarfixs.CURRENT_FIX_DURATION.max()


################ gets order for 2-3 target trials
tarordfix = tarfixs[tarfixs.ntargets.isin([2,3])]
tarordfix = tarordfix[tarordfix.refix==0]
order = np.ones(len(tarordfix), dtype=int)
for i in range(len(tarordfix)):
    if i > 0:
        if tarordfix.TRIAL_INDEX.iloc[i] == tarordfix.TRIAL_INDEX.iloc[i-1]:
            order[i] = 2
        if i > 1:
            if tarordfix.TRIAL_INDEX.iloc[i] == tarordfix.TRIAL_INDEX.iloc[i-1] & \
            tarordfix.TRIAL_INDEX.iloc[i] == tarordfix.TRIAL_INDEX.iloc[i-2]:
                order[i] = 3
tarordfix['Order'] = order 
#tarordfix.loc[:, 'fixation']   = 'target'
#tarordfix.to_csv('target_order_2_3.csv')


# remove the duplicates before rejoining
tarfixs = tarfixs[~tarfixs.fix_idx.isin(tarordfix.fix_idx.tolist())]
new = tarfixs.append(tarordfix)
new.loc[:, 'fixation']   = 'target'
new.to_csv('target_fixations.csv')
###############################################################################
# get distractor fixations
disfixs = fixs[fixs.CURRENT_FIX_RUN_INDEX == 1]
#disfixs = disfixs[disfixs.CURRENT_FIX_DURATION >= 100]
# only get fixations if they do not precede or follow a target fiaxtion
disfixs = disfixs[((disfixs.NEXT_FIX_INTEREST_AREA_LABEL == 'DISTRACTOR')|
                   (disfixs.NEXT_FIX_INTEREST_AREA_LABEL.isnull()))&
                  ((disfixs.PREVIOUS_FIX_INTEREST_AREA_LABEL == 'DISTRACTOR')|
                   (disfixs.PREVIOUS_FIX_INTEREST_AREA_LABEL.isnull()))]

                  
# get the distractors
disfixs = disfixs[disfixs.CURRENT_FIX_INTEREST_AREA_LABEL == 'DISTRACTOR']

disfixs = disfixs[disfixs.ntargets == 0]
disfixs = disfixs[disfixs.CURRENT_FIX_DURATION >= 120]
#disfixs = disfixs[disfixs.CURRENT_FIX_START > 1000]
#disfixs = disfixs[disfixs.CURRENT_FIX_START < disfixs.RT - 1500]

disfixs = disfixs[disfixs.CURRENT_FIX_DURATION < tar_max_dur]
disfixs.loc[:, 'fixation'] = 'distractor'
disfixs.to_csv('distractor_fixations.csv')

#%% check time between refixations
