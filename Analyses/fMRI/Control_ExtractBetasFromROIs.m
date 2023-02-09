% Parameter settings
clear; clc
%UserOptions.RoiType     = 'retinotopy'; %['retinotopy', 'loc', 'active cluster', 'deactive clusters'
UserOptions.RoiSpace    = 'native';

MaskName = {'Juelich_EVC', 'Juelich_VTC',...
    'Yeo_VAN', 'Yeo_DAN', 'Yeo_FPN', 'Yeo_DMN'};

UserOptions.ProjDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/fMRI';
%UserOptions.AnalysisDir = fullfile(UserOptions.ProjDir, 'Scripts');
%addpath(UserOptions.AnalysisDir)
UserOptions.TargetDir = fullfile(UserOptions.ProjDir, 'ROI_Data');

cope_list = [UserOptions.ProjDir '/COPEsForGoodRuns.mat'];
EVs = load(cope_list); EVs = EVs.EVs;

SJs = {'01' '04' '05' '07' '08' '09' '11' '13' '15' '16' '18' '19' '20'...
   '22' '25' '26' '29' '30' '31' '32' '33' '34' '35' '37' '38'};
%% Get data from ROIs from each condion in each block
for subj = 4:length(SJs) 
    SubEVs = EVs.(['sub' SJs{subj}]);
    for k = 1%:length(ROIs)
        [beta_hat, RoiSize] = Control_get_roi_data(SJs{subj}, SubEVs, UserOptions, MaskName{k});
        beta_hat(49) = [];          % discard no-response trials
        % put scr imgs to condition 41-48 
        scr_imgs = beta_hat(1,6:6:48); 
        beta_hat(6:6:48) = [];
        beta_hat(41:48) = scr_imgs;
        save(fullfile(UserOptions.TargetDir, ['sub' SJs{subj}], 'Betas', ['Control_' MaskName{k} '.mat']), 'beta_hat', 'RoiSize');
        clear  clear beta_hat real_beta RoiSize scr_beta Sw_hat u_hat z_data       
    end
    clear SubEVs
end
