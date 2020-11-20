# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 10:55:50 2016

@author: Joel
"""

import os

import pandas as pd
import numpy as np
from scipy import stats

# load fixations for all trials and subjects
fixs = pd.read_csv(
    '../data/fixation_report_all_trials_and_subjects.csv', na_values='.')

# add column to indicate whether it was a refixation
fixs['refix'] = 0
fixs.loc[fixs.CURRENT_FIX_REFIX_INTEREST_AREA != 0, 'refix'] = 1

# create a unique index for each fixation
fixs['fix_idx'] = fixs.index
fixs.index = fixs.timestamp.values

# remove practice trials and nonsense data
fixs = fixs.query('blockNum != 1 & TRIAL_INDEX <= 76')
fixs = fixs[fixs.ntargets.notnull()]

# list of target identity labels
target_lbls = ['TARGET_1','TARGET_2','TARGET_3']

# for each trial, determine how many unique targets have at least 1 fixation
fixs = fixs.set_index(['subsesh','TRIAL_INDEX']).sort_index()
fixs['ntargets_fixed'] = 0
for idx, df in fixs.groupby(level=['subsesh','TRIAL_INDEX']):
    if df.ntargets.iloc[0] == 0:
        fixs.loc[idx, 'ntargets_fixed'] = 0
    else:
        summary = df.CURRENT_FIX_INTEREST_AREA_LABEL.value_counts()
        summary = summary[summary.index.str.startswith('TARGET')]
        fixs.loc[idx, 'ntargets_fixed'] = summary.astype(bool).sum()
fixs.reset_index(inplace=True)  
  
# gather all target fixations whose location is within 1.5 degrees of the 
# center of the Interest Area 
tarfixs = fixs[fixs.CURRENT_FIX_NEAREST_INTEREST_AREA_LABEL.isin(target_lbls)]
                
# use only fixations which fall within 1.25 degrees of the target center, 
# which are not consecutive fixations in the same interest areas, which
# have a duration of at least 120 ms, which occur at least 1000 ms after the 
# trial started, which do not fall within 3000 ms of a button press
tarfixs = tarfixs.query('CURRENT_FIX_NEAREST_INTEREST_AREA_DISTANCE <= 1.25 \
                        & CURRENT_FIX_RUN_INDEX == 1 \
                        & CURRENT_FIX_DURATION >= 120 \
                        & CURRENT_FIX_START > 1000 \
                        & CURRENT_FIX_START  < RT - 3000')            

# additionally remove fixations whose duration is greater than the mean + 3SDs
# of all target fixations
tarfixs = tarfixs[(np.abs(stats.zscore(tarfixs.CURRENT_FIX_DURATION)) < 3)]
tar_max_dur = tarfixs.CURRENT_FIX_DURATION.max()

# get the order in which targets are fixated in trials with 2 or 3 targets
order_fixs = tarfixs[tarfixs.ntargets.isin([2.0,3.0])]
order_fixs = order_fixs[order_fixs.refix==0]
order = np.ones(len(order_fixs), dtype=int)
for i in range(len(order)):
    if i > 0:
        if order_fixs.TRIAL_INDEX.iloc[i] == order_fixs.TRIAL_INDEX.iloc[i-1]:
            order[i] = 2
        if i > 1:
            if (order_fixs.TRIAL_INDEX.iloc[i]==order_fixs.TRIAL_INDEX.iloc[i-1]
                & order_fixs.TRIAL_INDEX.iloc[i]==order_fixs.TRIAL_INDEX.iloc[i-2]):
                order[i] = 3
order_fixs['Order'] = order 
order_fixs.to_csv('../data/order_fixs_2_3_targets.csv')

# remove the duplicates before rejoining
tarfixs = tarfixs[~tarfixs.fix_idx.isin(order_fixs.fix_idx.tolist())]
tarfixs = tarfixs.append(order_fixs).sort_index()
tarfixs['fixation']   = 'target'
tarfixs.to_csv('../data/target_fixations.csv')

########################
# distractor fixations #
########################

# only look at distractor fixations if they do not precede or follow a target
# fiaxtion
disfixs = fixs[fixs.CURRENT_FIX_RUN_INDEX == 1]
disfixs = disfixs[((disfixs.NEXT_FIX_INTEREST_AREA_LABEL == 'DISTRACTOR')
                   | (disfixs.NEXT_FIX_INTEREST_AREA_LABEL.isnull())) 
                  & ((disfixs.PREVIOUS_FIX_INTEREST_AREA_LABEL == 'DISTRACTOR')
                   | (disfixs.PREVIOUS_FIX_INTEREST_AREA_LABEL.isnull()))]

# get the distractors
disfixs = disfixs.query('CURRENT_FIX_NEAREST_INTEREST_AREA_DISTANCE <= 1.25 \
                        & CURRENT_FIX_INTEREST_AREA_LABEL == "DISTRACTOR" \
                        & ntargets == 0 \
                        & CURRENT_FIX_DURATION >= 120 \
                        & CURRENT_FIX_START > 1000 \
                        & CURRENT_FIX_START  < RT - 3000')
                        
disfixs = disfixs[disfixs.CURRENT_FIX_DURATION < tar_max_dur]
disfixs['fixation'] = 'distractor'
disfixs.to_csv('../data/distractor_fixations.csv')
