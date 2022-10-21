#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 14:48:00 2022
- Compute the percentage of MEG-fMRI shared variance explained by a given model
  in a given brain network
- Carry out statistical comparisons between model in each brain network
- Plot the results (FIG 5 in MS)  
@author: wuy19
last update: 9/28/2022

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
FPCN = ['L MFG', 'R MFG', 'aPCC', 'R OFC']
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
# assign data to early vs late poststimulus period
bisect = np.zeros([df_all_poststimulus.shape[0],])
bisect[df_all_poststimulus['time'] <=1] =1
bisect[df_all_poststimulus['time'] >1] =2
df_all_poststimulus['Half'] =bisect

#separate df for each model 
df_RecImg = df_all_poststimulus[df_all_poststimulus['Model'] == 'Recognition']
df_TwoState = df_all_poststimulus[df_all_poststimulus['Model'] == 'TwoState']
#df for early and late period regardless models 
df_all_FirstHalf = df_all_poststimulus[df_all_poststimulus['Half'] ==1]
df_all_SecondHalf = df_all_poststimulus[df_all_poststimulus['Half'] ==2]
del bisect

# %% FIGURE 5 A-F
sns.set_style("white")
sns.set_context("paper")

for i, network in enumerate([EVC, VTC, SAL, DAN, FPCN, DMN]):
    data = df_all_poststimulus[df_all_poststimulus['ROI'].isin(network)] 
    mdata = data.groupby(['Model', 'time', 'Half'])['Explained variance'].agg(np.mean)
    mdata = mdata.to_frame()
    mdata.reset_index(inplace=True)
    
    fig, axes = plt.subplots(1,2, sharey= True, figsize=(4,2.2),
                         gridspec_kw={'width_ratios': [3,5]})    
    
    ax1= sns.pointplot(ax=axes[0], x='Half', y= 'Explained variance',
        hue='Model', estimator=np.median, alpha=0.5, ci=None,
        errwidth=1, capsize=0.1, marker = ".", data= mdata)
    
    #ax1.spines['left'].set_position(('outward', 5))
    ax1.yaxis.set_ticks_position('left')
    #ax1.spines['bottom'].set_position(('outward', 5))
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_bounds(0, 1)
    ax1.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax1.set_xticklabels(['Early', 'Late'], fontsize=8)
    ax1.set_yticklabels(['0', '25', '50', '75', '100'], fontsize=8)
    ax1.set_xlabel('Time', fontsize=8)
    ax1.set_ylabel('% of explained variance', fontsize=8)
    ax1.get_legend().remove() 
    
    ax2 = sns.violinplot(ax=axes[1], x='Half', y = 'Explained variance',
                    hue = 'Model', scale="count", inner="quartile",
                    data=mdata)
    ax2.yaxis.set_ticks_position('left')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_bounds(0, 1)
    ax2.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax2.set_xticklabels(['Early', 'Late'], fontsize=8)
    ax2.set_yticklabels(['0', '25', '50', '75', '100'], fontsize=8)
    ax2.set_xlabel('Time', fontsize=8)
    ax2.set_ylabel('')
    ax2.get_legend().remove()
    plt.show()
    plt.clf()
    del data, mdata, fig, axes, ax1, ax2
    
# FIGURE 5-SUPPLEMENT 1
for i, network in enumerate([EVC, VTC, SAL, DAN, FPCN, DMN]):
    fig, ax = plt.subplots(figsize=(2, 2))    
    ax = sns.pointplot(x='Half', y = 'Explained variance', hue = 'ROI',
                       estimator=np.median, alpha=0.5, ci=None, errwidth=1, 
                       capsize=0.1,linestyles='-', marker=None, 
        data = df_all_poststimulus[(df_all_poststimulus['ROI'].isin(network)) & \
                                  (df_all_poststimulus['Model']=='Recognition')])
    ax = sns.pointplot(x='Half', y = 'Explained variance',
        hue = 'ROI', estimator=np.median, alpha=0.5, ci=None,
        errwidth=1, capsize=0.1, linestyles='--', marker = None,
        data = df_all_poststimulus[(df_all_poststimulus['ROI'].isin(network)) & \
                                   (df_all_poststimulus['Model']=='TwoState')])                
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_bounds(0, 1)
    ax.set_yticks([-0.1, 0, 0.25, 0.5, 0.75, 1, 1.2])
    ax.set_xticklabels(['Early', 'Late'], fontsize=7)
    ax.set_yticklabels(['', '0', '25', '50', '75', '100', ''], fontsize=7)
    ax.set_xlabel('Time', fontsize=7)
    ax.set_ylabel('% of explained MEG-fMRI variance', fontsize=7)
    plt.title(printname[i], fontdict={'weight': 'bold', 'size': 7})
    plt.xticks(rotation = 0)
    ax.get_legend().remove()
    plt.show()
    plt.clf()

# %% DESCRIPTIVE STATISTICS & PAIR-WISE TESTS

Pvalues = []
def Q1(x):
    return np.percentile(x, 25)
def Q3(x):
    return np.percentile(x, 75)
    
for i, network in enumerate([EVC, VTC, SAL, DAN, FPCN, DMN]):
    data = df_all_poststimulus[df_all_poststimulus['ROI'].isin(network)] 
    mdata = data.groupby(['Model', 'time', 'Half'])['Explained variance'].agg(np.mean)
    mdata = mdata.to_frame()
    mdata.reset_index(inplace=True)  

    print('===================================================================')
    print('')
    print('')
    print(printname[i])
    print(mdata.groupby(['Model', 'Half'])['Explained variance'].agg([np.median, Q1, Q3]))
    print('')
# Mann-whitney test between the first and second half
    condition1 = mdata[(mdata['Model']=='Recognition') &  (mdata['Half']==1)]
    condition2 = mdata[(mdata['Model']=='Recognition') &  (mdata['Half']==2)]
    x = condition1['Explained variance']
    y = condition2['Explained variance']
    stats, pval = scipy.stats.mannwhitneyu(x,y,use_continuity=True, alternative='two-sided')
    print('Recognized Image 1 vs Recognized Image 2')
    print('statistics: ' + str(stats) + ', pval: ' + str(pval))
    print('')
    Pvalues.append(pval)
    del stats, pval, x, y, condition1, condition2

    condition1 = mdata[(mdata['Model']=='TwoState') &  (mdata['Half']==1)]
    condition2 = mdata[(mdata['Model']=='TwoState') &  (mdata['Half']==2)]
    x = condition1['Explained variance']
    y = condition2['Explained variance']
    stats, pval = scipy.stats.mannwhitneyu(x,y,use_continuity=True, alternative='two-sided')
    print('Two-State 1 vs Two-State 2')
    print('statistics: ' + str(stats) + ', pval: ' + str(pval))
    print('')
    Pvalues.append(pval)
    del stats, pval, x, y, condition1, condition2


#Wilcoxon test between models for each half respectively
    condition1 = mdata[(mdata['Model']=='Recognition') &  (mdata['Half']==1)]
    condition2 = mdata[(mdata['Model']=='TwoState') &  (mdata['Half']==1)]
    x = condition1['Explained variance']
    y = condition2['Explained variance']
    stats, pval = scipy.stats.wilcoxon(x,y, alternative='two-sided')
    print('Recognition 1 vs Two-State 1')
    print('statistics: ' + str(stats) + ', pval: ' + str(pval))
    print('')
    Pvalues.append(pval)
    del stats, pval, x, y, condition1, condition2

    condition1 = mdata[(mdata['Model']=='Recognition') &  (mdata['Half']==2)]
    condition2 = mdata[(mdata['Model']=='TwoState') &  (mdata['Half']==2)]
    x = condition1['Explained variance']
    y = condition2['Explained variance']
    stats, pval = scipy.stats.wilcoxon(x,y, alternative='two-sided')
    print('Recognized Image 2 vs Two-State 2')
    print('statistics: ' + str(stats) + ', pval: ' + str(pval))
    print('')
    Pvalues.append(pval)
    del stats, pval, x, y, condition1, condition2
    
rejected = multitest.multipletests(pvals=Pvalues, alpha=0.05, method='bonferroni')
#rejected = rejected[0]

# %% TWO-WAY MIXED-DESIGN ANOVA
import pingouin
for i, network in enumerate([EVC, VTC, SAL, DAN, FPCN, DMN]):
    
    data = df_all_poststimulus[df_all_poststimulus['ROI'].isin(network)] 
    mdata = data.groupby(['Model', 'time', 'Half'])['Explained variance'].agg(np.mean)
    mdata = mdata.to_frame()
    mdata.reset_index(inplace=True)  
    
    table = pd.pivot_table(data, values='Explained variance', index='Model',
                   columns='Half', aggfunc=np.median, margins=True)
    aov = pingouin.mixed_anova(data = mdata, subject = 'time', dv='Explained variance',
                     within= 'Model', between='Half', effsize='ng2', correction=True)
    print("")
    print('===================================================================')
    print(printname[i])
    print('===================================================================')
    print('')
    print(table)
    print('')
    print(aov)