# -*- coding: utf-8 -*-
"""
Created on Sun Oct 09 15:06:19 2016

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

#%% MAKE HDF STORES FOR THE DATA
sample_store_1000hz      = pd.HDFStore('samples_1000hz_r.h5')
sample_store_50hz        = pd.HDFStore('samples_50hz_r.h5')
range_store              = pd.HDFStore('ranges_r.h5')
trial_onset_store_BN_DS  = pd.HDFStore('trial_onset_BN_50hz_r.h5')
trial_offset_store_BN_DS = pd.HDFStore('trial_offset_BN_50hz_r.h5')

#%% MAIN PROCESSING LOOP

trial_info = ['ntargets','blockNum','CORRECT','subject','TRIAL_INDEX','RT']
start_msg, end_msg = 'STIM_ON', 'SEARCH_OVER'
messages = pd.read_table('message_report.txt')
messages['subsesh'] = messages.RECORDING_SESSION_LABEL
messages['subject'] = messages.RECORDING_SESSION_LABEL.str[:-2]
messages['session'] = messages.RECORDING_SESSION_LABEL.str[-1]
#start_msgs = messages[messages.CURRENT_MSG_TEXT=='STIM_ON']
#end_msgs   = messages[messages.CURRENT_MSG_TEXT=='SEARCH_OVER']
rp = pd.read_csv('reg_params.csv', index_col='subsesh')
rp = rp[rp.use==1]
rp.reset_index(inplace=True)
rp.set_index(['subject','session'], inplace=True)
rpmean = rp.mean()

for s in ss:
    print 'Processing: ' + str(s)
    sub, sesh = int(s[1:-2]), int(s[-1:])
    ascfile = datadir + s + os.sep + s + '.asc'
    s_starts   = events_from_message_report(messages, s, start_msg)
    s_ends     = events_from_message_report(messages, s, end_msg)
    samps, events = load_eyelink_dataset(ascfile)
    samps         = interp_eyelink_blinks(samps, events, find_recovery=True, interp_fields='pup_l')
    samps         = interp_zeros(samps, interp_fields=['pup_l']) 
    
    if s == 'p20_2':
        samps['c_pup_l'] = (samps.pup_l - (samps.x_l * rpmean.x_l_b) - (samps.y_l * rpmean.y_l_b).values)
    else:
        samps['c_pup_l'] = (samps.pup_l - (samps.x_l * rp.loc[sub,sesh].x_l_b) - (samps.y_l * rp.loc[sub,sesh].y_l_b).values)

    
    samps         = butterworth_series(samps,filt_order=3,cutoff_freq=0.008, fields=['pup_l','c_pup_l'])
    samps50hz     = samps[::20].copy()
    ranges        = extract_events(samps, 
                                   s_starts,
                                   offset=0,
                                   duration=10000,
                                   units='samples',
                                   borrow_attributes=trial_info)
    baselines     = extract_events(samps, 
                                   s_starts,
                                   offset=-1000,
                                   duration=1000,
                                   units='samples',
                                   borrow_attributes=trial_info).mean(level=0)
    trialon       = extract_events(samps50hz, 
                                   s_starts,
                                   offset=-150,
                                   duration=300,
                                   units='samples',
                                   borrow_attributes=trial_info)                               
    trialoff      = extract_events(samps50hz, 
                                   s_ends,
                                   offset=-150,
                                   duration=300,
                                   units='samples',
                                   borrow_attributes=trial_info)
    
    ranges['pup_l_change']   = (((ranges.pup_l - baselines.pup_l) / baselines.pup_l) * 100).values # % signal change from baseline
    trialon['pup_l_change']  = (((trialon.pup_l - baselines.pup_l) / baselines.pup_l) * 100).values # % signal change from baseline
    trialoff['pup_l_change'] = (((trialoff.pup_l - baselines.pup_l) / baselines.pup_l) * 100).values # % signal change from baseline
    ranges['c_pup_l_change']   = (((ranges.c_pup_l - baselines.c_pup_l) / baselines.c_pup_l) * 100).values # % signal change from baseline
    trialon['c_pup_l_change']  = (((trialon.c_pup_l - baselines.c_pup_l) / baselines.c_pup_l) * 100).values # % signal change from baseline
    trialoff['c_pup_l_change'] = (((trialoff.c_pup_l - baselines.c_pup_l) / baselines.c_pup_l) * 100).values # % signal change from baseline

    print_info(events, samps, ranges)    
    
    sample_store_50hz.put(s, samps50hz)
    sample_store_1000hz.put(s, samps)           
    range_store.put(s, ranges)  
    trial_onset_store_BN_DS.put(s, trialon)
    trial_offset_store_BN_DS.put(s, trialoff)

sample_store_50hz.close()        
sample_store_1000hz.close()   
range_store.close()
trial_onset_store_BN_DS.close()
trial_offset_store_BN_DS.close()

#%% SUBTRACTIVE BASELINE CORRECTION
range_store              = pd.HDFStore('ranges_r.h5')
trial_onset_store_BN_DS  = pd.HDFStore('trial_onset_BN_50hz_r.h5')
trial_offset_store_BN_DS = pd.HDFStore('trial_offset_BN_50hz_r.h5')

for s in ss:
    ranges   = range_store[s]
    trialon  = trial_onset_store_BN_DS[s]
    trialoff = trial_offset_store_BN_DS[s]
    trialon.reset_index(inplace=True)
    baselines = trialon.loc[((trialon.onset>=100)&(trialon.onset<150))]
    baselines.set_index(['event','onset'], inplace=True)
    baselines = baselines.mean(level=0)
    ranges['c_pup_l_sb']   = (ranges.c_pup_l - baselines.c_pup_l).values # subtractive baseline  
    trialon['c_pup_l_sb']  = (trialon.c_pup_l - baselines.c_pup_l).values # subtractive baseline  
    trialoff['c_pup_l_sb'] = (trialoff.c_pup_l - baselines.c_pup_l).values # subtractive baseline       
    range_store.put(s, ranges)  
    trial_onset_store_BN_DS.put(s, trialon)
    trial_offset_store_BN_DS.put(s, trialoff)
range_store.close()
trial_onset_store_BN_DS.close()
trial_offset_store_BN_DS.close()
#%% LINEAR WARPING OF PUPIL DATA DURING SEARCH
range_store.open()
normalized_ranges = pd.HDFStore('normalized_ranges_r.h5')
trial_info = ['interpolated','ntargets','blockNum','CORRECT','subject','TRIAL_INDEX','RT']

subs = list(set([val[:-2] for val in ss]))

all_subs_df = pd.DataFrame()    
for s in subs:
    ranges_1 = range_store[s + '_1']
    ranges_2 = range_store[s + '_2']
    normal_ranges_1 = normalise_time(ranges_1, trial_info)
    normal_ranges_2 = normalise_time(ranges_2, trial_info)
    normal_ranges_1.reset_index(inplace=True)
    normal_ranges_2.reset_index(inplace=True)
    normal_ranges_2.event = normal_ranges_2.event + 76
    #normal_ranges_2.set_index(['event','onset'])
    norm_ranges = normal_ranges_1.append(normal_ranges_2)
    normalized_ranges.put(s, norm_ranges)
    all_subs_df = all_subs_df.append(norm_ranges)
    print 'Done: ' + str(s)

#all_subs_df['session'] = all_subs_df['subject'].str[-1:] 
all_subs_df['subject'] = all_subs_df['subject'].str[1:]  
all_subs_df.to_csv('all_subs_norm_ranges_r_sb.csv')

range_store.close()
normalized_ranges.close()

#%% get non-self-term-data
range_store.open()
subs = list(set([val[:-2] for val in ss]))
func = {'interpolated':lambda x: float(x.sum()) /float(x.count()),'c_pup_l_sb':'mean'}
df = pd.DataFrame()
for s in subs:
    ranges_1 = range_store[s + '_1']
    ranges_2 = range_store[s + '_2']
    ranges_1.reset_index(inplace=True)
    ranges_2.reset_index(inplace=True)
    ranges_2.event = ranges_2.event + 76
    ranges = ranges_1.append(ranges_2)
    ranges = ranges.loc[ranges.RT=='10000']
    ranges = ranges.loc[ranges.CORRECT=='1']
    ranges.interpolated = ranges.interpolated.convert_objects(convert_numeric=True)
    df = df.append(ranges.groupby(['subject','ntargets','onset'], as_index=False).aggregate(func))

df = df[::20]
df.to_csv(datadir + "search_non_self_term_ranges_sb.csv")

#%% SORT OUT ONSET AND OFFSET DATA
trial_onset_store_BN_DS.open()
trial_offset_store_BN_DS.open()

subs = list(set([val[:-2] for val in ss]))

onpup_df  = pd.DataFrame()
offpup_df = pd.DataFrame()
for s in subs:
    on1  = trial_onset_store_BN_DS[s + '_1']
    on2  = trial_onset_store_BN_DS[s + '_2']
    off1 = trial_offset_store_BN_DS[s + '_1']
    off2 = trial_offset_store_BN_DS[s + '_2']
    for df in [on1,on2,off1,off2]:
        df.reset_index(inplace=True)
    for df in [on2,off2]:
        df.event = df.event + 76
    onpupdat  = on1.append(on2)
    offpupdat = off1.append(off2)
    onpupdat['subject'] = onpupdat['subject']
    offpupdat['subject'] = offpupdat['subject']
    onpup_df  = onpup_df.append(onpupdat)
    offpup_df = offpup_df.append(offpupdat)

onpup_df.to_csv('all_subs_onpup_r_sb.csv')
offpup_df.to_csv('all_subs_offpup_r_sb.csv')

trial_onset_store_BN_DS.close()
trial_offset_store_BN_DS.close()

#%% GET FIXATION-ALIGNED RESPONSES (TAKES ABOUT 12-13 HOURS)
fix_aligned_store = pd.HDFStore('fix_aligned_50hz_r_100msb.h5')
sample_store      = pd.HDFStore('samples_1000hz_r.h5')
#tarfixs_1_2 = pd.read_csv('target_fixations_1_2.csv', index_col='timestamp')
tarfixs     = pd.read_csv('target_fixations.csv', index_col='timestamp')
disfixs     = pd.read_csv('distractor_fixations.csv', index_col='timestamp')
#tarfixs_3   = pd.read_csv('target_fixations_3.csv', index_col='timestamp')
#tarfixs_3   = pd.read_csv('target_order_2_3.csv', index_col='timestamp')


trial_info    = ['subject',
                 'session',
                 'refix',
                 'Order',
                 'CORRECT',
                 'CURRENT_FIX_DURATION',
                 'CURRENT_FIX_START',
                 'TRIAL_INDEX',
                 'ntargets',
                 'ntargets_fixed',
                 'blockNum',
                 'fixation',
                 'fix_idx',
                 'RT']

for s in ss:
    print 'Processing ' + str(s)
    samps     = sample_store[s]
    #s_tar_1_2 = tarfixs_1_2[tarfixs_1_2.subsesh==s]
    s_tar = tarfixs[tarfixs.subsesh==s]
    s_dis = disfixs[disfixs.subsesh==s]
    #s_tar_3   = tarfixs_3[tarfixs_3.subsesh==s]
    print 'Getting fixaligned for ' + str(len(s_tar)) + ' target fixations'
    tar_ranges = get_fix_aligned(samps, s_tar, trial_info)
    print 'Getting fixaligned for ' + str(len(s_dis)) + ' distractor fixations'
    dis_ranges = get_fix_aligned(samps, s_dis, trial_info)
    #print 'Getting fixaligned for ' + str(len(s_tar_3)) + ' Order fixations'
    #tar_ranges_3   = get_fix_aligned(samps, s_tar_3, trial_info)
    #df = pd.concat([tar_ranges_1_2, dis_ranges, tar_ranges_3])
    df = pd.concat([tar_ranges, dis_ranges])

    df = df[::20]
    #allonsets = [range(125) for i in range(len(df)/125)]
    #allonsets = [val for sublist in allonsets for val in sublist]
    df['onset'] = df['onset'] / 20
    df.reset_index(inplace=True)
    fix_aligned_store.put(s, df)
    print 'Done: ' + str(s)


fix_aligned_store.close()
sample_store.close()
#%% COMPUTE SUBJECT AVERAGES

fix_aligned_store = pd.HDFStore('fix_aligned_50hz_r_100msb.h5')

subs = set([val[:-2] for val in ss])

grand_df_012 = pd.DataFrame()
grand_df_3   = pd.DataFrame()

for s in subs: 
    dats_1 = fix_aligned_store[s+'_1']
    dats_2 = fix_aligned_store[s+'_2']
    dat_1_2 = pd.concat([dats_1, dats_2])
    
    nt012     = dat_1_2[dat_1_2.Order.isnull()]
    nt012.fixation.loc[nt012.refix==1] = nt012.fixation.loc[nt012.refix==1] + '_refix'
    nt012     = nt012[nt012.CORRECT==1]
#    nt012     = nt012[nt012.CURRENT_FIX_START > 1000]
#    nt012     = nt012[nt012.CURRENT_FIX_DURATION >= 150]
#    nt012 = nt012[nt012.CURRENT_FIX_START <= nt012.RT - 1500]

    #nt012     = nt012[nt012.CURRENT_FIX_START < nt012.RT]
    nt012_agg = nt012.groupby(['subject','fixation','onset'], as_index=False)['pup_l_change','c_pup_l_change'].mean()
    grand_df_012 = grand_df_012.append(nt012_agg)
    
    orderfixs = dat_1_2[~dat_1_2.Order.isnull()]
    orderfixs = orderfixs[orderfixs.CORRECT==1]
#    orderfixs = orderfixs[orderfixs.CURRENT_FIX_DURATION >= 150]
#    orderfixs = orderfixs[orderfixs.CURRENT_FIX_START <= orderfixs.RT - 1500]
#
#    orderfixs = orderfixs[orderfixs.CURRENT_FIX_START >= 1000]
    orderfixs_agg = orderfixs.groupby(['subject','Order','onset'], as_index=False)['pup_l_change','c_pup_l_change'].mean()
    grand_df_3 = grand_df_3.append(orderfixs_agg)

grand_df_012.to_csv('fixaligned_012_r_100msb.csv')
grand_df_3.to_csv('fixaligned_3_r_100msb.csv')

#%%
fix_aligned_store.open()
all_fixs = pd.DataFrame()
for s in ss:
    dat = fix_aligned_store[s]
    #dat = dat[dat.onset >=5]
    #new = dat.groupby(['subject','fix_idx'],as_index=False)['pup_l','CURRENT_FIX_DURATION'].mean()
    all_fixs = all_fixs.append(dat)
    
all_fixs.to_csv('all_fixs_50hz_r_100msb.csv')    

#%% GET EACH SUBJECTS PERCENTAGE OF INTERPOLATED DATA
#range_store              = pd.HDFStore('ranges_no_r.h5')
range_store              = pd.HDFStore('ranges_r.h5')

subs = set([val[:-2] for val in ss])

pc_int_d = {}
for sub in subs:
    sesh1 = range_store[sub + '_1']
    sesh2 = range_store[sub + '_2']
    int_s1 = sesh1.interpolated.value_counts(normalize=True)[1]
    int_s2 = sesh2.interpolated.value_counts(normalize=True)[1]
    pc_int = (int_s1 + int_s2) / 2
    pc_int_d[sub] = pc_int

df = pd.DataFrame(pc_int_d, index=[0])
range_store.close()

#%%

td  = pd.read_table('trialreport.csv', na_values = '.', sep=',')
td.subject = td.subject.str.strip('p')
td.subject = pd.to_numeric(td.subject)
range_store              = pd.HDFStore('ranges_r.h5')

for key in range_store.keys():
    sub  = int(key[2:-2])
    sesh = int(key[-1])
    rs = range_store.get(key)
    rs.interpolated = pd.to_numeric(rs.interpolated)
    pc_interp = rs.groupby(level=[0])['interpolated'].mean()*100
    td.loc[(td.subject==sub)&(td.session==sesh), 'pct_interp'] = pc_interp.values
    
td.to_csv('new_trialreport.csv')

#%%

td = pd.read_csv('new_trialreport.csv')
aggregations = {'RT':'mean',
                'CORRECT':'mean',
                'selfterm':'mean',
                'BLINK_COUNT':'mean',
                'FIXATION_COUNT':'mean',
                'SACCADE_COUNT':'mean',
                'AVERAGE_BLINK_DURATION':'mean',
                'AVERAGE_FIXATION_DURATION':'mean',
                'PUPIL_SIZE_MEAN':'mean',
                'VISITED_INTEREST_AREA_COUNT':'mean'}
                
ag_td = td.groupby(['subject','ntargets']).agg(aggregations)

#%%
import matplotlib.pyplot as plt
import seaborn as sn
df = pd.DataFrame()
for key in range_store.keys():
    dat = range_store[key]
    dat = dat.mean(level=1)
    dat['ss'] = key
    dat.reset_index(inplace=True)
    df = df.append(dat)