% Computing averaged 1-Pearson's r dissimilarity distancess between
% BOLD-signals to presented images, for each subject and ROI separately.
% We first compute the mean beta parameter estimate across all trials
% corresponding to the same image, resulting in 40 mean beta parameter 
% vectors. Then, we compute the 1-Pearson's r between all possible mean
% beta vector pairs.

TrialType = 'real'; %['real', 'scr', 'all'];

DistanceMeasure  = 'correlation';
DataDir          = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/fMRI/ROI_Data'; % this is where the beta parameter patterns are stored
TargetDir= 'Unnormalized_Distances';

if strcmp(TrialType, 'all') ==1
   conds = 1:48;
elseif strcmp(TrialType, 'real') ==1
    conds = 1:40;
elseif strcmp(TrialType, 'scr') ==1
    conds = 41:48;
end

nConditions = length(conds);
nDistances  = (nConditions*nConditions-nConditions)/2;

SJs = {'01' '04' '05' '07' '08' '09' '11' '13' '15' '16' '18' '19' '20'...
    '22' '25' '26' '29' '30' '31' '32' '33' '34' '35' '37' '38'};

ROIs = {'V1', 'V2', 'V3', 'V4', 'loc_face', 'loc_animal', 'loc_house', 'loc_object', ...
    'active_Brainstem', 'active_BG_Thalamus', 'active_IPS_L', ...
    'active_IPS_R', 'active_aPCC', 'active_aInsula_L', 'active_aInsula_R', ...
    'active_IFJ_R', 'active_MFG_L', 'active_MFG_R', 'active_OFC_R',...
    'deactive_AG_L', 'deactive_AG_R', 'deactive_HC_L', 'deactive_HC_R', ...
    'deactive_mPFC', 'deactive_PCC', 'deactive_SFG_L', 'deactive_SFG_R', ...
    'deactive_STG_L', 'deactive_STG_R'};
%% ************************************************************************
for r = 1:length(ROIs)
    for subj = 1%:length(SJs)
        % load subject-specific pattern vectors for a  given ROI
        SubjData = fullfile(DataDir,['sub' SJs{subj}], 'Betas', [ROIs{r} '.mat']);
        SubjData = load(SubjData);
        SubjData = SubjData.beta_hat;
        
        % check if data exist. If not, continue with the next subject
        if iscell(SubjData) ~= 1
            AvgDistances.(ROIs{r}) = nan;
            save(fullfile(DataDir, ['sub' S/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/fMRI/ROI_Data/sub01/Unnormalized_Distances/ControlAnalysisJs{subj}], TargetDir,...
                'allAvg1-Pearson.mat'), 'AvgDistances')
            clear subjData
            disp('                                                             ')
            disp('*************************************************************')
            disp(['Analysis on sub' SJs{subj} ' ' ROIs{r} ' done'])
            disp('*************************************************************')
            disp('                                                             ')
            continue
        else
        end
        
        % if data exist, compute the mean beta estimate for each stimulus
        nVoxels = size(SubjData{1,1},2);
        AvgBeta_hat = zeros(nConditions, nVoxels);
        for n = conds
            AvgBeta_hat(n,:) = nanmean(SubjData{1,n},1);
        end
        % compute 1-Pearson's r distances between all  pairs and save them
        AvgDistances.(ROIs{r}) = pdist(AvgBeta_hat, DistanceMeasure);
        
        save(fullfile(DataDir, ['sub' SJs{subj}], TargetDir,...
            'realAvg1-Pearson.mat'), 'AvgDistances')
        clear AvgBeta_hat SubjData nVoxels
        disp('                                                             ')
        disp('*************************************************************')
        disp(['Analysis on sub' SJs{subj} ' ' ROIs{r} ' done'])
        disp('*************************************************************')
        disp('                                                             ')
    end
end
