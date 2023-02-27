#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 10:59:06 2023

@author: wuy19
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

RootDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/eLife_Wu2022/'
imgs_dict = mat73.loadmat(RootDir + 'Figures/Data/n_trials_per_image.mat')
nTrials =  imgs_dict['n_trials']
del imgs_dict



face_r = nTrials[:,np.arange(0,5,1)]
house_r = nTrials[:,np.arange(5,10,1)]
object_r = nTrials[:,np.arange(10,15,1)]
animal_r = nTrials[:,np.arange(15,20,1)]
face_u = nTrials[:,np.arange(20,25,1)]
house_u = nTrials[:,np.arange(25,30,1)]
object_u = nTrials[:,np.arange(30,35,1)]
animal_u = nTrials[:,np.arange(35,40,1)]
place_holder = np.zeros([24,1])

data = np.concatenate((face_r, place_holder, house_r, place_holder,
                       object_r, place_holder, animal_r, place_holder,
                       face_u, place_holder, house_u, place_holder,
                       object_u, place_holder, animal_u), axis=1)
mean_nTrials = np.mean(data,axis=0)
sem_nTrials = np.std(data, ddof=1, axis=0) / np.sqrt(np.size(nTrials, axis=0))


green =  '#2ca02c' #44 160 44
red = '#d62728'
colors = np.tile('#2ca02c', (24,)).tolist() + np.tile('#d62728', (23,)).tolist() 
condition_names = ['face R', 'house R', 'object R', 'animal R',
              'face U', 'house U', 'object U', 'animal U']


fig, ax = plt.subplots(figsize=(8,2))
ax.bar(np.arange(0,len(mean_nTrials),1), mean_nTrials, yerr=sem_nTrials, width =0.8,
       capsize=2, color=colors)
ax.set_yticks(np.arange(0,20,4))
ax.set_xticks(np.arange(2,48,6))
ax.set_xticklabels(condition_names, fontsize=9)
ax.set_ylabel('Mean trial number per image', fontsize=9)
ax.spines['left'].set_position(('outward', 5))
ax.yaxis.set_ticks_position('left')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.vlines(5,0,16, colors='k', linestyles = 'dotted')
ax.vlines(11,0,16, colors='k', linestyles = 'dotted')
ax.vlines(17,0,16, colors='k', linestyles = 'dotted')
ax.vlines(23,0,16, colors='k', linestyles = 'dotted')
ax.vlines(23,0,16, colors='k', linestyles = 'dotted')
ax.vlines(29,0,16, colors='k', linestyles = 'dotted')
ax.vlines(35,0,16, colors='k', linestyles = 'dotted')
ax.vlines(41,0,16, colors='k', linestyles = 'dotted')
plt.xticks(rotation = 0)
plt.tight_layout()
plt.savefig('TrialsPerExemplars.svg')
#plt.show