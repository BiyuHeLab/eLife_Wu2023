#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon, Jan 30, 2023

- serve as additional control analysis for the percentage of MEG-fMRI shared
  variance explained by a given model in a given brain network
- Carry out statistical comparisons between model in each brain network
- divide data over 2 sec into 10, 20, or 40 time bins

@author: wuy19
last update: 2/16/2023

"""

import mat73
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import scipy
from statsmodels.stats import multitest
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'


ROIs = ['V1', 'V2', 'V3', 'loc_face', 'loc_animal', 'loc_house', 'loc_object',
    'active_IPS_L', 'active_IPS_R', 'active_MFG_L', 'active_MFG_R', 
    'active_IFJ_R', 'active_aPCC',
    'active_aInsula_L', 'active_aInsula_R', 'active_OFC_R',
    'deactive_AG_L', 'deactive_AG_R',
    'deactive_mPFC', 'deactive_PCC', 'deactive_SFG_L', 'deactive_SFG_R',
    'deactive_STG_L', 'deactive_STG_R']

PlotLabels = ['V1', 'V2', 'V3', 'face', 'animal', 'house', 'object',
    'L IPS','R IPS', 'L MFG', 'R MFG', 'R IFJ', 'aPCC',
    'L aInsula', 'R aInsula', 'R OFC',
    'L AG', 'R AG', 'mPFC', 'PCC', 'L SFG', 'R SFG',
    'L STG', 'R STG']

#Assign ROIs to brain network (Yeo et al., 2012)
EVC = ['V1', 'V2', 'V3']
VTC = ['animal', 'face', 'house', 'object']
SAL = ['L aInsula', 'R aInsula'] 
DAN = ['R IFJ', 'L IPS', 'R IPS']    
FPCN = ['L MFG', 'R MFG', 'R OFC','aPCC']
DMN = ['L AG', 'R AG', 'PCC', 'mPFC',  'L SFG', 'R SFG']
printname = ['EVC', 'VTC', 'SAL', 'DAN', 'FPCN', 'DMN']

# %% DATA PREPARATION


comm_dict = mat73.loadmat('./Data/Fusion_Commonality.mat')
fusion_dict = mat73.loadmat('./Data/Fusion_MEG-MRI_Correlation.mat');
comm_Rec = comm_dict['C']['Recognition']
comm_TwoState = comm_dict['C']['TwoState']
fusion = fusion_dict['AvgRho']
del comm_dict, fusion_dict

Stats = mat73.loadmat('./Data/Fusion_RecognitionModel_ClusterInference.mat')
Stats_Rec= Stats['ClusterInference']['Recognition']
Stats = mat73.loadmat('./Data/Fusion_TwoStateModel_ClusterInference.mat')
Stats_TwoState = Stats['ClusterInference']['TwoState']
del Stats

# Merge information from different .mat files to a single pd.DataFrame file
df = pd.DataFrame()
for model in ["Recognition", "TwoState"]:

    if model == "Recognition":
       data = comm_Rec
       sig = Stats_Rec
    elif model =="TwoState":
       data = comm_TwoState
       sig = Stats_TwoState
      
    for i, r in enumerate(ROIs):        
        df_ROI = pd.DataFrame(data= np.arange(-0.5,2.0025, 0.0025), columns = ['time'])
        df_ROI['ROI'] = PlotLabels[i]
        df_ROI['Model'] = model
        df_ROI['Commonality Coefficient'] = data[r]
        df_ROI['Fusion'] = np.power(fusion[r],2)
        tmp = data[r]
        tmp[tmp <= 0] = 0
        df_ROI['Explained variance'] = tmp / np.power(fusion[r],2)  
        df_ROI['Significance'] = sig[r]['SigTimePoint_SumPos']
        df = pd.concat([df, df_ROI])
        del df_ROI, tmp
    del data, sig
del Stats_Rec, Stats_TwoState, comm_Rec, comm_TwoState, i, r

# select poststimulus only
df_all_poststimulus = df[(df['time'] >=0)]
n_split = 10
time_idx = np.kron(np.arange(1,n_split+1,1), np.ones([1, int(800/n_split)]))
time_idx = np.append(time_idx, n_split)
time_idx = np.tile(time_idx, (48))

df_all_poststimulus['Split'] =time_idx
del time_idx
del fusion, model
# %% FIGURE 5 A-F
sns.set_style("white")
sns.set_context("paper")

for i, network in enumerate([EVC, VTC, SAL, DAN, FPCN, DMN]):
    data = df_all_poststimulus[df_all_poststimulus['ROI'].isin(network)] 
    mdata = data.groupby(['Model', 'time', 'Split'])['Explained variance'].agg(np.mean)
    mdata = mdata.to_frame()
    mdata.reset_index(inplace=True)
    
    #Wilcoxon test between models for each half respectively
    Pvals = np.zeros(n_split,)
    Pvals[:]=np.nan
    for t, split in enumerate(np.arange(1,n_split+1,1)):
        condition1 = mdata[(mdata['Model']=='Recognition') &  (mdata['Split']==split)]
        condition2 = mdata[(mdata['Model']=='TwoState') &  (mdata['Split']==split)]
        x = condition1['Explained variance']
        y = condition2['Explained variance']
        stats, pval = scipy.stats.wilcoxon(x,y, alternative='two-sided')
        Pvals[t] = pval
        rejected = multitest.multipletests(pvals=Pvals, alpha=0.05, method='fdr_bh')   
        del condition1, condition2, x,y, stats, pval
    
    fig, ax = plt.subplots(figsize=(9,1.5))    
    ax= sns.pointplot(x='Split', y= 'Explained variance',
    hue='Model', estimator=np.median, alpha=0.5, errorbar= ('pi', 50),
    errwidth=None, capsize=0.1, marker = None, data= mdata)
        
    ax.spines['left'].set_position(('outward', 5))
    ax.yaxis.set_ticks_position('left')
    ax.spines['bottom'].set_position(('outward', 5))
    ax.xaxis.set_ticks_position('bottom')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_bounds(0, 1)
    ax.set_xticks(np.linspace(-0.5,n_split-0.5, n_split+1))
    xticklabels = np.linspace(0,2,n_split+1)
    xticklabels = np.around(xticklabels,2)
    ax.set_xticklabels(xticklabels.astype(str), fontsize=9)
    ax.set_xlabel('Time (s) relative to stimulus oneset', fontsize=9)
    
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1, 1.25])
    ax.set_yticklabels(['0', '25', '50', '75', '100', ''], fontsize=9)  
    ax.set_ylabel('% of explained variance', fontsize=9)
    ax.get_legend().remove()
    ax.set_title(printname[i], fontweight='bold', fontsize=12)
    
    for k in range(0,n_split):
        if (rejected[1][k] < 0.05) and (rejected[1][k] > 0.01):
            ax.annotate('*', (k-0.2, 1.1), fontsize=9)
        elif (rejected[1][k] < 0.01) and (rejected[1][k] > 0.001):
            ax.annotate('**', (k-0.2, 1.1), fontsize=9)
        elif (rejected[1][k] < 0.001):
            ax.annotate('***', (k-0.2, 1.1), fontsize=9)
    plt.savefig('ModelComparison_' + str(n_split) + '_' + printname[i] + '.svg', dpi=600, bbox_inches='tight')        
    del ax, fig, rejected, k, data, mdata, Pvals,t 
