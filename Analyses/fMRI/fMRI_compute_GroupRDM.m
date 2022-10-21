clc;clear;

ProjDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/';
DataDir = fullfile(ProjDir, 'fMRI/ROI_Data');
addpath(fullfile(ProjDir, 'fMRI/Scripts'))
nConditions = 40;

SJs = {'01', '04', '05', '07', '08', '09', '11',  '13', '15', '16', '18', '19',...
    '20', '22', '25', '26', '29', '30', '31', '32', '33', '34', '35', '37', '38'};

measurement = 'realAvg1-Pearson';
Normalization = 'Unnormalized';
AvgDir = [DataDir '/Grand_Average/' Normalization '_Data'];

ROIs =  {'V1', 'V2', 'V3', 'V4', 'loc_face', 'loc_animal', 'loc_house', 'loc_object',...
    'active_Brainstem', 'active_BG_Thalamus', 'active_IPS_L', ...
    'active_IPS_R', 'active_aPCC', 'active_aInsula_L', 'active_aInsula_R',...
    'active_IFJ_R', 'active_MFG_L', 'active_MFG_R', 'active_OFC_R',...
    'deactive_AG_L', 'deactive_AG_R', 'deactive_HC_L', 'deactive_HC_R',...
    'deactive_mPFC', 'deactive_PCC', 'deactive_SFG_L', 'deactive_SFG_R',...
     'deactive_STG_L', 'deactive_STG_R'};


nDistances = (nConditions * nConditions - nConditions)/2;
%%

for r = 1:length(ROIs)
    tmp = nan(length(SJs), nDistances);
    for sj = 1:length(SJs)
        SubData = fullfile(DataDir, ['sub' SJs{sj}], 'Unnormalized_Distances', [measurement '.mat']);
        SubData = load(SubData);
        roinames = fieldnames(SubData.AvgDistances);
        if ~ismember(ROIs{r}, roinames)
            continue
        else           
            roi_data = SubData.AvgDistances.(ROIs{r});
      
            if size(roi_data,2)~=1
                tmp(sj,:) = roi_data;
            else
                tmp(sj,:) = nan(1,nDistances);
            end
            clear roi_data SubData
        end
    end
    AvgData.(ROIs{r}) = nanmean(tmp,1);
    clear tmp
end
save(fullfile(AvgDir, [measurement '.mat']), 'AvgData')