# eLife_Wu2022

Analysis scripts and codes for reproducing the figures for the manuscript Spatiotemporal neural dynamics of objeccts recognition under uncertainty in humans**

Analysis was originally performed on a machine running Rehat Linux v7.8, Matlab R2017a, and Spyder 5.0.5 (Python 3.8).

## Contents

### **Figures/**
data and code to reproduce the figures

**Figures/data/** contains all analysis outputs required to plot the figure 1,2 4, and 5 in the manuscript. 

* bhv_df.pkl contains behavioral responses in the fMRI experiment.
* percent_seen_real contains MEG subjects' recognition rates in recognized real images trials.
* correct_seen.p contains MEG subjects' categorization accuracies in recognized real images trials.
* correct_unseen.p contains MEG subjects' categorization accuracies in unrecognized real images trials.

* MEG_decoding_Recognition.mat contains MEG recognition outcome decoding result (Figure 2A)
* MEG_decoding_Recognition_shuffled.mat contains permutation statistics necessary for MEG recognition decoding cluster inference.  
* MEG_decoding_Recognition_stats.mat contains cluster-level statistics for MEG recognition outcome decoding.    

* MEG_decoding_Object_R.mat contains MEG obejct category decoding result for real recognized trials (Figure 2B, green)
* MEG_decoding_Object_U.mat contains MEG obejct category decoding result for real unrecognized trials (Figure 2B, red)
* MEG_decoding_Object_shuffled.mat contains permutation statistics necessary for MEG object category decoding cluster inference.  
* MEG_decoding_Object_stats.mat contains cluster-level statistics for MEG object category decoding.

* MEG_decoding_Localizer_.mat conatains obejct category decoding result for the MEG localizer block (Figure 2-supplement 1)
* MEG_decoding_Localizer_shuffled.mat contains permutation statistics necessary for MEG object category localizer decoding cluster inference. 

* MEG_avgRDM_100Hz.mat contains group-average MEG RDMs for time points corresponding to Figure 2C and D     

* Fusion_Commonality.mat contains commonality coefficients for all ROIs included in the analysis (Figure 4 and 5).  
* Fusion_MEG-MRI Correlation.mat contains the bivariate correlations (Spearman's rho) between MEG and fMRI RDMs 
* Fusion_RecognitionModel_ClusterInference.mat contains cluster-level statistics for commonality coefficients corresponding to the recognition model  
* Fusion_TwoStateModel_ClusterInference.mat contains cluster-level statistics for commonality coefficients corresponding to the recognition model 

### **Analyses/**
Codes to output the presented results

**Analyses/MEG/**

* HLTP.py: paramteres and helper functions for MEG preprocessing
* MEG_preprocessing.py: detrend, deamean, band-pass filter, subsample, epoch ICA cleaned MEG data
* MEG_MNE2fieldtrip.m: transforming data from MNE to fieldtrip format
* MEG_compute_distances.m: calculate 1-r distances between pairs of images
* MEG_compute_GroupRDM.m: compute the group average RDM for visualization and commonality analysis
* MEG_prepareDataForDecoding.m: reorganize the data for decoding analysis 
* MEG_decoding_RecognitionOutcomes.m: perform recognition outcome decoding
* MEG_decoding_RecognitionOutcomes_ClusterInference.m: perform cluster inference for recognition outcome decoding
* MEG_decoding_object.m: perform object category decoding in the main task
* MEG_decoding_object_ClusterInfernce.m: perform cluster inference for object category decoding in the main task
* import_mne_epochs.m: subfunction for MEG_MNE2fieldtrip.m
* organize_ERF.m: subfunction for MEG_MNE2fieldtrip.m

**Anlayses/fMRI/**
* fMRI_preprocessing_RunFirstLevelFeat: perfom GLM at subject level using FSL
* fMRI_preprocessing_StimIdentityGLM_template.fsf: GLM parameter-specifiction
* fMRI_ExtractBetasFromROI.m: extract and store the beta estimates from each ROI
* fMRI_compute_distances.m: compute 1-r distances between pairs of images
* fMRI_conpute_GroupRDM.m: compute group-average RDMs for commonality analysis
* fMRI_get_roi_data.m: subfunction for fMRI_ExtractBetaFromROI.m
* load_nifi.m: subfunction for fMRI_ExtractBetaFromROI.m

**Analyses/Fusion**
* Fusion_CommonalityAnalysis.m: main script to run commonality analysis and generate the permuation samples
* Create_ModelVectors.m: subfunction for creating model RDMs
* compute_commonality.m: subfunction for executing commonality analysis
* correlation.m: subfunction for computing semipartial correlations (@author: Martin Hebart (2018, eLife)
