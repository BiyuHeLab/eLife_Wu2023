% Parameter settings
clear; clc
UserOptions.RoiType     = 'retinotopy'; %['retinotopy', 'loc', 'active cluster', 'deactive clusters'
UserOptions.RoiSpace    = 'native';

if strcmp(UserOptions.RoiType, 'retinotopy')==1
    MaskName = {'V1', 'V2', 'V3', 'V4'};
    ROIs =     {'V1', 'V2', 'V3', 'V4'};
    
elseif strcmp(UserOptions.RoiType, 'loc')==1
    MaskName = {'animal_2mm_maskedVTC_noEVC_Z2.3_500_voxels_highres.nii',...
        'face_2mm_maskedVTC_noEVC_Z2.3_500_voxels_highres.nii',...
        'house_2mm_maskedVTC_noEVC_Z2.3_500_voxels_highres.nii',...
        'object_2mm_maskedVTC_noEVC_Z2.3_500_voxels_highres.nii'};
    ROIs =      {'loc_animal', 'loc_face', 'loc_house', 'loc_object'};
    
elseif strcmp(UserOptions.RoiType, 'active clusters')==1
    MaskName = {'Brainstem_LR.nii', 'Caudate_Putamen_Thalamus_LR.nii', ...
        'IPS_LOC_Precuneus_L.nii', 'IPS_R.nii', 'Posterior_Cingulate_LR.nii',...
        'Anterior_Insula_L.nii', 'Anterior_Insula_R.nii', 'IFJ_R.nii',...
        'Frontal_Gyri_L.nii', 'Frontal_Gyri_R.nii', 'OFC_R.nii'};    
    ROIs =  {'Brainstem', 'BG_Thalamus', 'IPS_L', ...
        'IPS_R', 'aPCC', 'aInsula_L', 'aInsula_R', 'IFJ_R',...
        'MFG_L', 'MFG_R', 'OFC_R'};

elseif strcmp(UserOptions.RoiType, 'deactive clusters')==1
    MaskName = {'Angular_Gyrus_L.nii', 'Angular_Gyrus_R.nii',...
        'Hippocampus_L.nii', 'Hippocampus_R.nii', 'Medial_PFC_LR.nii', ...
        'Precuneus_LR.nii', 'Superior_Frontal_Gyrus_L.nii', ...
        'Superior_Frontal_Gyrus_R.nii', ...
        'Temporal_Gyri_L.nii', 'Temporal_Gyri_R.nii'};
    ROIs = {'AG_L', 'AG_R', 'HC_L', 'HC_R', 'mPFC', 'PCC', 'SFG_L', 'SFG_R', ...
        'STG_L', 'STG_R'};
end

UserOptions.ProjDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/fMRI';
UserOptions.AnalysisDir = fullfile(UserOptions.ProjDir, 'Scripts');
addpath(UserOptions.AnalysisDir)
UserOptions.TargetDir = fullfile(UserOptions.ProjDir, 'ROI_Data');

cope_list = [UserOptions.ProjDir '/COPEsForGoodRuns.mat'];
EVs = load(cope_list); EVs = EVs.EVs;

SJs = {'01' '04' '05' '07' '08' '09' '11' '13' '15' '16' '18' '19' '20'...
   '22' '25' '26' '29' '30' '31' '32' '33' '34' '35' '37' '38'};
%% Get data from ROIs from each condion in each block
for subj = 1:length(SJs) 
    SubEVs = EVs.(['sub' SJs{subj}]);
    for k = 1%:length(ROIs)
        [beta_hat, RoiSize] = fMRI_get_roi_data(SJs{subj}, SubEVs, UserOptions, MaskName{k});
         beta_hat(49) = [];          % discard no-response trials
        % put scr imgs to condition 41-48 
        scr_imgs = beta_hat(1,6:6:48); 
        beta_hat(6:6:48) = [];
        beta_hat(41:48) = scr_imgs;
        save(fullfile(UserOptions.TargetDir, ['sub' SJs{subj}], 'Betas', ['active_' ROIs{k} '.mat']), 'beta_hat', 'RoiSize');
        clear  clear beta_hat real_beta RoiSize scr_beta Sw_hat u_hat z_data       
    end
    clear SubEVs
end
