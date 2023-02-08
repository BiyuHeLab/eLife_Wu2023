#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 29 14:25:39 2021
Combine multiple ROIs to a network masks
1. vmPFC:  
2. DMN
3. Salience
4. Visual  

@author: wuy19
"""
ProjDir= '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/fMRI/'
import sys
sys.path.insert(0, ProjDir)
import os 
os.chdir(ProjDir)

import HLTP
from nipype.interfaces import fsl 
import numpy as np
import nibabel
import nilearn
from nilearn import plotting, image, surface
from nilearn import datasets
from nilearn import plotting, surface, image, datasets, input_data, reporting
from nilearn.input_data import NiftiMasker

mask_fname ='/usr/local/fsl/data/standard/MNI152_T1_2mm_brain_mask.nii.gz'
template_fname ='/usr/local/fsl/data/standard/MNI152_T1_2mm_brain.nii.gz'
mask_img = nibabel.load(mask_fname)
template = nibabel.load(template_fname)
# %% create V1, V2, V3v, V3d, and V4 mask from Juelich atlas
mask_folder = ProjDir + 'ROI_masks/MPM/'
roi_idx = 113
roi_imgs = []

# select ROIs being combined
for i, hemi in enumerate(['l', 'r']):
    roi_img = mask_folder + 'JulichBrain_MPMAtlas_' + hemi + '_N10.nii.gz'
    
    roi_2mm = image.resample_img(roi_img, target_affine=mask_img._affine,
                            target_shape=[91,109,91])
    roi_2mm._dataobj[roi_2mm._dataobj !=roi_idx]=0
    roi_2mm._dataobj[roi_2mm._dataobj ==roi_idx]=1
    roi_imgs.append(roi_2mm)
    del roi_img, roi_2mm
    
bilateral_img = nilearn.image.math_img('img1 + img2',
                                   img1=roi_imgs[0],
                                   img2=roi_imgs[1])
bilateral_img._data[bilateral_img._data !=0] =1
roi_size = len(bilateral_img._data[bilateral_img._data ==1])

plotting.plot_roi(bilateral_img, bg_img=template, colorbar = False,
                      draw_cross=False, display_mode='ortho')
bilateral_img.to_filename(mask_folder + 'V4_bl.nii.gz')
del bilateral_img, roi_imgs, roi_idx, roi_size


rois =  ['V1', 'V2', 'V3d', 'V3v', 'V4']
roi_imgs = []

for i in range(len(rois)):
    roi_img = nibabel.load(mask_folder + rois[i] + '_bl.nii.gz')
    roi_imgs.append(roi_img)
    del roi_img

EVC_img = nilearn.image.math_img('img1 + img2 + img3 + img4 + img5',
                                   img1=roi_imgs[0],
                                   img2=roi_imgs[1],
                                   img3=roi_imgs[2],
                                   img4=roi_imgs[3],
                                   img5=roi_imgs[4])
EVC_img._data[EVC_img._data !=0] = 1
EVC_img.to_filename(ProjDir + 'ROI_masks/Masks_ControlAnalysis/EVC_standard.nii.gz')
plotting.plot_roi(EVC_img, bg_img=template, colorbar = False,
                      draw_cross=False, display_mode='ortho')
del roi_imgs, rois
# %% 1) Create mask in ventral visual stream based on HArvard Oxford atlas (MNI space)
#   2.) merge these mask to a single ROI mask
mask_folder = ProjDir + 'ROI_masks/Harvard_ROIs/'
roi_idx = {'LOCs': 22, 'LOCi': 23, 'PPH': 35, 'TFCp': 38, 'TOF': 39, 'OFG': 40}
#roi_imgs = []
roi_size= []

# select ROIs being combined
for key, value in roi_idx.items():
    print(key, value)
    roi_img = nibabel.load(mask_folder + 'HarvardOxford-cort-maxprob-thr25-2mm.nii')
    roi_map = roi_img.get_fdata()
    roi_map[roi_map !=value]=0
    roi_map[roi_map ==value]=1
    bilateral_img = nilearn.image.new_img_like(roi_img, roi_map)     
    #roi_imgs.append(roi_2mm)
    #del roi_img, roi_2mm
    bilateral_img.dataobj[bilateral_img.dataobj !=0] =1
    roi_size.append(len(bilateral_img._data[bilateral_img._data ==1]))
    plotting.plot_roi(bilateral_img, bg_img=template, colorbar = False,
                      draw_cross=False, display_mode='ortho')
    bilateral_img.to_filename(mask_folder + key + '_bl.nii.gz')
    del bilateral_img, roi_map, roi_img
    
#combine higher visual areas, but no LOCs     
rois =  ['LOCi', 'PPH', 'TFCp', 'TOF', 'OFG']
roi_imgs = []

for i in range(len(rois)):
    roi_img = nibabel.load(mask_folder + rois[i] + '_bl.nii.gz')
    roi_imgs.append(roi_img)
    del roi_img

VTC_img = nilearn.image.math_img('img1 + img2 + img3 + img4 + img5',
                                   img1=roi_imgs[0],
                                   img2=roi_imgs[1],
                                   img3=roi_imgs[2],
                                   img4=roi_imgs[3],
                                   img5=roi_imgs[4])
plotting.plot_roi(VTC_img, bg_img=template, colorbar = False,
                      draw_cross=False, display_mode='ortho')
                      
# %% find overlapping voxels between EVC and VTC parcels and remove them from
# the VTC parcel                      
EVC_img._data[EVC_img._data !=0] = 3

plotting.plot_roi(EVC_img, bg_img=template, colorbar = False,
                      draw_cross=False, display_mode='ortho')

merged_img = nilearn.image.math_img('img1 + img2',
                                   img1=EVC_img,
                                   img2=VTC_img)
                             
merged_img._data[merged_img._data !=1] = 0
len(merged_img._data[merged_img._data ==1])
plotting.plot_roi(merged_img, bg_img=template, colorbar = False,
                      draw_cross=False, display_mode='ortho')       
merged_img.to_filename(mask_folder + 'VTC_bl.nii.gz')
 
# %% Split the Yeo_native_mask to 7 masks in each subject's native space
roi_folder = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/fMRI/ROI_masks/'

import HLTP
roi_dict = {'Visual':1, 'Somatosensory':2 , 'DAN': 3, 'VAN': 4, 'Limbic': 5, 'FPN': 6, 'DMN': 7}

DataDir = '/isilon/LFMI/VMdrive/data/HLTP_fMRI'
RegPath = '/proc_data/func/block1/StimIdentityGLM.feat/reg'



for subj in HLTP.subjects[1:]:
    SubDir = roi_folder + 'sub' + str(subj).zfill(2) + '/'
   
    
    for key, value in roi_dict.items():
        Yeo_img = nibabel.load(SubDir + 'Yeo7_native.nii.gz')        
        roi_map = Yeo_img.get_fdata()
        roi_map[roi_map != value]=0
        roi_img = nilearn.image.new_img_like(Yeo_img, roi_map)  
        roi_img.to_filename(SubDir + 'Yeo_' + key + '_native.nii.gz')
        del roi_map, roi_img, Yeo_img
    del SubDir

# %% SPlit the MNI Yeo mask to seven separate masks in MNI space.
Yeo_img = nibabel.load(roi_folder + 'Yeo_7_networks_2mm.nii.gz')
Yeo_img._dataobj = np.squeeze(Yeo_img._dataobj, axis=3)

for key, value in roi_dict.items():        
    roi_map = Yeo_img.get_fdata()
    roi_map[roi_map != value]=0
    roi_img = nilearn.image.new_img_like(mask_img, roi_map)  
    roi_img.to_filename(roi_folder + 'Yeo_' + key + '_standard.nii.gz')
    del roi_map, roi_img

# %% Transfoorm the Julich-EVC and VTC masks from MNI to native space, for each subject 

roi_folder = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/fMRI/ROI_masks/Masks_ControlAnalysis/'
DataDir = '/isilon/LFMI/VMdrive/data/HLTP_fMRI'
RegPath = '/proc_data/func/block1/StimIdentityGLM.feat/reg'

rois = ['EVC', 'VTC']

for subj in HLTP.subjects:
    ant_vol = DataDir + '/sub' + str(subj).zfill(2) + \
        '/proc_data/anat/divt1pd_brain_2mm_mask.nii'
    RegFile = DataDir + '/sub' + str(subj).zfill(2) + RegPath + \
        '/standard2highres.mat'  
    SubMaskDir = ProjDir + '/ROI_masks/sub' + str(subj).zfill(2) + '/'
    
    for label in rois:
        flt = fsl.FLIRT()
        flt.inputs.in_file = roi_folder + label + '_standard.nii.gz'
        flt.inputs.reference = ant_vol
        flt.inputs.out_file = SubMaskDir + 'Juelich_' + label + '_native.nii.gz'
        flt.inputs.output_type = 'NIFTI_GZ'
        flt.inputs.interp = 'nearestneighbour'
        flt.inputs.apply_xfm = True
        flt.inputs.in_matrix_file = RegFile
        flt.cmdline
        flt.run()
        del flt
    del ant_vol, RegFile, SubMaskDir
 