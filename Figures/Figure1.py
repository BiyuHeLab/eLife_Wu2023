#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:44:48 2022

Statistics and plots for recognition and categorization behavior in MEG and 
fMRI experiment (FIG 1C and D)  
@author: wuy19
LAst Update: 10/20/2022
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
##############################################################################
# STATISTICS ON BEHAVIOR
############################################################################## 

#subset data by recognition report
MEG_seen = pd.read_pickle('./Data/percent_seen_real.p')
MEG_unseen = pd.read_pickle('./Data/percent_seen_real.p')   

fMRI_bhv = pd.read_pickle('./Data/bhv_df.pkl')   
proportion_R = fMRI_bhv.groupby(['real', 'subject'])['R'].mean()
fMRI_seen = 100. * proportion_R.xs(1, level='real').values
del proportion_R    
    
df_MEG = pd.DataFrame(data=MEG_seen, columns= ['performance'])
df_MEG['modality'] = 'MEG'
df_fMRI = pd.DataFrame(data=fMRI_seen, columns= ['performance'])
df_fMRI['modality'] = 'fMRI'
df_report = pd.concat([df_MEG, df_fMRI])


# Statistics for Figure 1C
rate = df_MEG['performance'].mean()
sem = (df_MEG['performance'].std(ddof=1)) / np.sqrt(df_MEG['performance'].shape[0])
w, p = stats.wilcoxon(x=df_MEG['performance']-50, alternative="two-sided")
print("")
print("MEG RECOGNITION RATE : " + str(np.round(rate,1)) + "% SEM: " + str(np.round(sem,1)))
print("Wilcoxon signed rank test against chance: W:" + str(w), ", p-val: " + str(p))  
print("")


rate = df_fMRI['performance'].mean()
sem = (df_fMRI['performance'].std(ddof=1)) / np.sqrt(df_fMRI['performance'].shape[0])
w, p = stats.wilcoxon(x=df_fMRI['performance']-50, alternative="two-sided")
print("")
print("fMRI RECOGNITION RATE : " + str(np.round(rate,1)) + "% SEM: " + str(np.round(sem,1)))
print("Wilcoxon signed rank test chance: W:" + str(w), ", p-val: " + str(p))
print("")


    
# subset data by categorization behavior 
MEG_correct_R = pd.read_pickle('./Data/correct_seen.p')
MEG_correct_U = pd.read_pickle('./Data/correct_unseen.p')
df_MEG_R = pd.DataFrame(data=MEG_correct_R, columns= ['performance'])
df_MEG_R['modality'] = 'MEG'
df_MEG_R['report'] = 'R'
df_MEG_U = pd.DataFrame(data=MEG_correct_U, columns= ['performance'])
df_MEG_U['modality'] = 'MEG'
df_MEG_U['report'] = 'U'

      
correct_pd_group = fMRI_bhv.groupby(['recognition', 'real', 'subject'
                                   ])['correct'].mean()
fMRI_correct_R = 100. * correct_pd_group.xs((1, 1), 
                                    level=('recognition','real')).values
fMRI_correct_U = 100. * correct_pd_group.xs((-1, 1), 
                                    level=('recognition','real')).values
df_fMRI_R = pd.DataFrame(data=fMRI_correct_R, columns= ['performance'])
df_fMRI_R['modality'] = 'fMRI'
df_fMRI_R['report'] = 'R'
df_fMRI_U = pd.DataFrame(data=fMRI_correct_U, columns= ['performance'])
df_fMRI_U['modality'] = 'fMRI'
df_fMRI_U['report'] = 'U'
df_cat = pd.concat([df_MEG_R, df_MEG_U, df_fMRI_R, df_fMRI_U])


# Statistics in Figure 1D: MEG
rate = df_MEG_R['performance'].mean()
sem = (df_MEG_R['performance'].std(ddof=1)) / np.sqrt(df_MEG_R['performance'].shape[0]) 
w, p = stats.wilcoxon(x=df_MEG_R['performance']-25, alternative="greater")
print("")
print("MEG categorization accuracy for recognized real images  : " + str(np.round(rate,1)) + "% SEM: " + str(np.round(sem,1)))
print("Wilcoxon signed rank test chance: W:" + str(w), ", p-val: " + str(p))
print("")

rate= df_MEG_U['performance'].mean()
sem = (df_MEG_U['performance'].std(ddof=1)) / np.sqrt(df_MEG_U['performance'].shape[0])  
w, p = stats.wilcoxon(x=df_MEG_U['performance']-25, alternative="greater")
print("")
print("MEG Recognition accurcay for unrecognized real images  : " + str(np.round(rate,1)) + "% SEM: " + str(np.round(sem,1)))
print("Wilcoxon signed rank test chance: W:" + str(w), ", p-val: " + str(p)) 
print("")

w, p = stats.wilcoxon(x=df_MEG_R['performance'], y=df_MEG_U['performance'], alternative="greater")
print("")
print("Difference in categorization accuracy between MEG R and MEG U ")
print("Wilcoxon: W:" + str(w), ", p-val: " + str(p)) 
print("")


# Statistics in Figure 1D: fMRI
rate = df_fMRI_R['performance'].mean()
sem = (df_fMRI_R['performance'].std(ddof=1)) / np.sqrt(df_fMRI_R['performance'].shape[0]) 
w, p = stats.wilcoxon(x=df_fMRI_R['performance']-25, alternative="greater")
print("")
print("fMRI categorization accuracy for recognized real images  :  : " + str(np.round(rate,1)) + "% SEM: " + str(np.round(sem,1)))
print("Wilcoxon signed rank test chance: W:" + str(w), ", p-val: " + str(p))
print("")

rate = df_fMRI_U['performance'].mean()
sem = (df_fMRI_U['performance'].std(ddof=1)) / np.sqrt(df_fMRI_U['performance'].shape[0])
print("")  
w, p = stats.wilcoxon(x=df_fMRI_U['performance']-25, alternative="greater")
print("fMRI Recognition accurcay for unrcognized real images  : " + str(np.round(rate,1)) + "% SEM: " + str(np.round(sem,1)))
print("Wilcoxon signed rank test chance: W:" + str(w), ", p-val: " + str(p))
print("")

w, p = stats.wilcoxon(x=df_fMRI_R['performance'], y=df_fMRI_U['performance'], alternative="greater")
print("")
print("Difference in categorization accuracy between fMRI R and fMRI U ")
print("Wilcoxon: W:" + str(w), ", p-val: " + str(p)) 
print("")
del df_MEG_R, df_MEG_U, df_fMRI_R, df_fMRI_U, df_MEG, df_fMRI, rate, sem, w, p

# %% PLOT THE DATA

sns.set_style("white")
sns.set_context("paper")
fig, axes = plt.subplots(1,2, sharey= True, figsize=(3.5,2),
                         gridspec_kw={'width_ratios': [2,2.8]})   

ax1 = sns.violinplot(ax=axes[0], x='modality', y = 'performance', data=df_report,
                    palette='pastel', inner='quartile',
                    scale = 'count', linewidth=0.5) # cut = 0)
ax1 = sns.swarmplot(ax=axes[0], x ='modality', y ='performance', data = df_report,
                   edgecolor='gray', palette='dark', size = 3)

ax1.spines['left'].set_position(('outward', 5))
ax1.yaxis.set_ticks_position('left')

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['left'].set_bounds(0, 100)
ax1.set_yticks([0, 25, 50, 75, 100])
ax1.set_yticklabels(['0', '25', '50', '75', '100'], fontsize=8)
ax1.set_xlabel('', fontsize=7)
ax1.set_ylabel('% Recognized', fontsize=8)
ax1.axhline(y=50, xmin=0, xmax=1, c="k", linestyle="--" )

ax2 = sns.violinplot(ax=axes[1], x='modality', y = 'performance', hue='report', data=df_cat,
                    palette='pastel', inner='quartile',
                    scale = 'count', linewidth=0.5) # cut = 0)
ax2 = sns.swarmplot(ax=axes[1], x ='modality', y ='performance', hue='report', data = df_cat,
                   edgecolor='gray', palette='dark', dodge=True, size=3)

ax2.spines['left'].set_position(('outward', 5))
ax2.yaxis.set_ticks_position('left')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['left'].set_bounds(0, 100)
ax2.set_yticks([0, 25, 50, 75, 100])

ax2.set_yticklabels(['0', '25', '50', '75', '100'], fontsize=8)
ax2.set_xlabel('', fontsize=7)
ax2.set_ylabel('% Correct categorization', fontsize=8)
ax2.axhline(y=25, xmin=0, xmax=1, c="k", linestyle="--" )  
ax2.get_legend().remove()
plt.tight_layout()


