%==========================================================================
% Main script to compute commonality coefficient of a single model, MRI and
% MEG data, using the actual and shuffled data
% REQUIRES:   'CommonalityAnalysis_SingleModel.m'
%             'Create_ModelVectors.m'
%             'Correlate.m'               
%
% Author: @Yuan-hao Wu
%Latest Update: 9/25/2022 @Yuan-hao Wu
%==========================================================================
clear; clc

addpath('./helper functions')
fMRI_data = 'path to fMRI RDM';
MEG_data = 'path to MEG RDM';

ROIs = fieldnames(fMRI_data.AvgData);
ModelNames = {'Recognition', 'TwoState'};
    
times = -0.5:0.0025:2;
%% Commonality analysis for each model and ROI, respetively
for m = 1:length(ModelNames)
    modelRDM = Create_ModelVectors(ModelNames{m});
    megRDM = MEG_data.AvgData'; 
    
    for r = 1:length(ROIs)
        fmriRDM = fMRI_data.AvgData.(ROIs{r})';
        commModel  = compute_commonality(megRDM,fmriRDM, modelRDM);
        C.(ModelNames{m}).(ROIs{r})= commModel;
        clear commModel fmriRDM
    end
    clear modelRDM
end
save('path_output_file', 'C')
%% Shuffling fMRI-RDM image order and compute shared variance between the model, MEG, and shuffled fMRI

nPerms = 5000;
C_perm = nan(nPerms, length(times),length(ROIs));
modelRDM = Create_ModelVectors(ModelNames{1});
megRDM = MEG_data.AvgData';
parfor r = 1:length(ROIs)
    tic
    for k = 1:nPerms
        fmriRDM = fMRI_data.AvgData.(ROIs{r})';
        fmriRDM = fmriRDM(randperm(length(fmriRDM)));
        C_perm(k,:,r) = CommonalityAnalysis_SingleModel(megRDM,fmriRDM, modelRDM);
    end
    
    elapsed = toc;
    disp('                                                        ')
    disp('==========================================================================')
    disp([ModelNames{1} '-' ROIs{r} ' done. Elpased time: ' num2str(elapsed) ' sec.'])
    disp('==========================================================================')
    disp('                                                        ')
    
end
clear modelRDM megRDM k m r elapsed
save('path_output_file', 'C_perm')

