# eLife_Wu2022

Analysis scripts and codes for reproducing the figures for the manuscript Spatiotemporal neural dynamics of objeccts recognition under uncertainty in humans**

Analysis was originally performed on a machine running Rehat Linux v7.8, Matlab R2017a, and Spyder 5.0.5 (Python 3.8).

## Contents

### **Figures/**
data and code to reproduce figures

**Figures/data/** contains all analysis outputs required to plot the figure 1,2 4, and 5 in the manuscript. 

*
*bhv_df.pkl* contains behavioral responses in the fMRI experiment.
*percent_seen_real* contains MEG subjects' recognition rates in recognized real images trials.
*correct_seen.p* contains MEG subjects' categorization accuracies in recognized real images trials.
*correct_unseen.p* contains MEG subjects' categorization accuracies in unrecognized real images trials.

*MEG_decoding_Recognition.mat* contains MEG recognition outcome decoding result (Figure 2A)
*MEG_decoding_Recognition_shuffled.mat* contains permutation statistics necessary for MEG recognition decoding cluster inference.  
*MEG_decoding_Recognition_stats.mat* contains cluster-level statistics for MEG recognition outcome decoding.    

*MEG_decoding_Object_R.mat* contains MEG obejct category decoding result for real recognized trials (Figure 2B, green)
*MEG_decoding_Object_U.mat* contains MEG obejct category decoding result for real unrecognized trials (Figure 2B, red)
*MEG_decoding_Object_shuffled.mat* contains permutation statistics necessary for MEG object category decoding cluster inference.  
*MEG_decoding_Object_stats.mat* contains cluster-level statistics for MEG object category decoding.

*MEG_decoding_Localizer_.mat* conatains obejct category decoding result for the MEG localizer block (Figure 2-supplement 1)
*MEG_decoding_Localizer_shuffled.mat* contains permutation statistics necessary for MEG object category localizer decoding cluster inference. 

*MEG_avgRDM_100Hz.mat* contains group-average MEG RDMs for time points corresponding to Figure 2C and D     

*Fusion_Commonality.mat* contains commonality coefficients for all ROIs included in the analysis (Figure 4 and 5).  
*Fusion_MEG-MRI Correlation.mat* contains the bivariate correlations (Spearman's rho) between MEG and fMRI RDMs 
*Fusion_RecognitionModel_ClusterInference.mat* contains cluster-level statistics for commonality coefficients corresponding to the recognition model  
*Fusion_TwoStateModel_ClusterInference.mat* contains cluster-level statistics for commonality coefficients corresponding to the recognition model 







analysis/
code to analyze the data

Ana
