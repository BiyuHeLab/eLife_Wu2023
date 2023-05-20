#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 18:59:29 2020

"""
import sys
sys.path.insert(0, '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/eLife_Wu2022/Analyses')
HomeDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/proc_data/'
import os
import mne
from scipy import signal
import HLTP


for SUB_ID in HLTP.subjects:
    SubDir = HomeDir + SUB_ID 
    tmin, tmax = -3, 3
    block_epochs = [];
    subj_raw_dir = HLTP.MEG_raw_dir + '/' + SUB_ID + '/'
    #fdir = HLTP.MEG_pro_dir + '/' + SUB_ID
    _, _, n_blocks, _ = HLTP.get_experimental_details(subj_raw_dir) 


    for b in range(1, n_blocks + 1):
        fname = SubDir + '/thresh' + str(b).zfill(2) + '_stage2_raw.fif'
        raw = mne.io.read_raw_fif(fname, preload=True)
    
    #*****************************************************************************
    # Remove linear trend in the raw data
    #*****************************************************************************
        raw.info['bads'] = HLTP.bads[SUB_ID]  
        picks = mne.pick_types(raw.info, meg = True, ref_meg = False,
                               exclude = 'bads')
        raw.apply_function(signal.detrend, picks=picks, dtype=None) # n_jobs=24)
             
    #*****************************************************************************
    # Apply a low-pass filter of 35 Hz to the data
    #*****************************************************************************          
  
        raw.filter(l_freq=None, h_freq=35, picks=picks, 
                   filter_length='auto', phase = HLTP.filter_phase, 
                   method = HLTP.filter_method, h_trans_bandwidth = 'auto'), #n_jobs=24)
          
    #*****************************************************************************
    # Epoch the continuous data
    #*****************************************************************************  
        events = mne.find_events(raw, stim_channel=HLTP.stim_channel, 
                                 mask = HLTP.event_id['Stimulus'], 
                                 mask_type = 'and')  
    # Correct for photo-diode delay:
        events[:, 0] = events[:, 0] + HLTP.PD_delay
    
        decim = 3            
        epochs = mne.Epochs(raw, events, {'Stimulus': 1}, tmin=tmin, tmax=tmax, 
                            proj=False, baseline=None, picks=picks, preload=True, detrend=0, 
                            verbose=False, decim= decim)
        block_epochs.append(epochs) 
        del epochs, events, fname, raw
    
    for b in range(0, n_blocks):
        block_epochs[b].info['dev_head_t'] = block_epochs[0].info['dev_head_t']
   
    
    all_epochs = mne.concatenate_epochs(block_epochs)
    del block_epochs
    all_epochs.event_id = HLTP.event_dict
    conditions = HLTP.get_conditions(subj_raw_dir, n_blocks)
    all_epochs.events[:,2] = conditions
   
    bad_trials = HLTP.load(HLTP.MEG_pro_dir + '/' + SUB_ID + '/bad_trials.p')
    if len(bad_trials) < 40:
       all_epochs.drop(bad_trials)
   
    all_epochs.save(os.path.join(SubDir, 'HLTP_0-35Hz_400Hz-epo.fif'))
    del all_epochs, bad_trials, conditions, n_blocks, picks, subj_raw_dir, SubDir, b
