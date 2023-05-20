% =========================================================================
% Compute correlation between the object category model and
% between-exemplar distances for recognized and unrecognized trials,
% respectively.
% The category model RDM for a given rocgnition outcome has a size of
% 20x20 cells (corresponding to 4 object categrories x 5 exemplars).

% input data:
% Subject's MEG RDV in n_timepoints x n_conditionpairs x n_subjects format
% n_conditionpairs = 780 which corresponds to a 40x40 RDM

% outputs:
% MEG_RSA_CategoryModel.mat (Spearman's rho between MEG and model RDM)
% MEG_RSA_CategoryModel_stats.mat (Stats correspond to the cluster-based
% permutation test)

% Last update: 3/7/2023 @Yuan-hao Wu
% =========================================================================

clear; clc;
addpath('../../Figures/helper functions')
OutputDir = 'PATH FOR THE OUTPUT FILES';
SubjData = load('ADD PATH FOR THE REQUIRED INPUT DATA');

times = -0.5:0.01:2;
conditions = {'seen', 'unseen'};
n_subjects = size(SubjData,3);
%% %Create 20 x20 Category Model RDM
CategoryRDM =ones(20,20);
for i = [1 6 11 16]
    CategoryRDM(i:i+4, i:i+4) = 0;
end
for i = 1:20
    CategoryRDM(i,i) = 0;
end
ModelVector = squareform(CategoryRDM)';
%% Compute Spearman's rhos between  time-varying MEG data - category model
% for each recognition outcome, respectively.
% output variable: rho (struct with two fields)

SubRDMs.seen = zeros(190, length(times), n_subjects);
SubRDMs.unseen = zeros(190, length(times), n_subjects);

for sub = 1:n_subjects
   
    for t_idx = 1:length(times)
        tmp = squareform(SubjData(t_idx,:,sub));
        SubRDMs.seen(:, t_idx, sub) = squareform(tmp(1:20,1:20));
        SubRDMs.unseen(:, t_idx, sub) = squareform(tmp(21:40,21:40));        
        clear tmp 
    end
    
    [rho1, ~] = corr(SubRDMs.seen(:,:,sub), ModelVector, 'type', 'Spearman', 'rows', 'complete', 'tail', 'right');
    [rho2, ~] = corr(SubRDMs.unseen(:,:,sub), ModelVector, 'type', 'Spearman', 'rows', 'complete', 'tail', 'right');   
    rho.seen(sub,:)= rho1;
    rho.unseen(sub,:)= rho2;
    
    clear rho1 rho2
end
save(fullfile(OutputDir, 'MEG_RSA_CategoryModel.mat'), 'rho');
%% statistical test against zero correlation at each time points using
% one-taield Wilcoxon test
% output variable: Wilcoxon (struct with 2 fields)

for i_time = 1:length(times)
    
    [p1,h,stats] = signrank(rho.seen(:,i_time), 0, 'tail','right');
    Wilcoxon.seen.p(1,i_time) = p1; Wilcoxon.seen.h(1,i_time) = h;
    Wilcoxon.seen.zval(1,i_time) = stats.zval; Wilcoxon.seen.signedrank(1,i_time) = stats.signedrank;
    clear p h stats
    
    [p1,h,stats] = signrank(rho.unseen(:,i_time), 0, 'tail','right');
    Wilcoxon.unseen.p(1,i_time) = p1; Wilcoxon.unseen.h(1,i_time) = h;
    Wilcoxon.unseen.zval(1,i_time) = stats.zval; Wilcoxon.unseen.signedrank(1,i_time) = stats.signedrank;
    clear p h stats
end
%% Generate empirical null distribution for the model effect & do cluster inference
n_perms = 5000;
[PermRho, PermStats] = MEG_RSA_CategoryModel_Permutations(ModelVector, SubRDMs, n_perms);

CDT = 0.05;
ClusterStats = 'SumPos';  % ['SumPos', 'ClusterSize]

for c = 1:length(conditions)
    clusters_orig = find_temporal_clusters(Wilcoxon.(conditions{c}).signedrank , Wilcoxon.(conditions{c}).p, CDT);
    ClusterInference.(conditions{c}).clusters_orig = clusters_orig;
    clear clusters_orig
    
    for i = 1:n_perms
        cluster_shuffled = find_temporal_clusters(PermStats.(conditions{c}).signedrank(i,:),PermStats.(conditions{c}).p(i,:), CDT);
        
        ClusterInference.(conditions{c}).maxClusterSizes(1,i) = cluster_shuffled.maxSize;
        ClusterInference.(conditions{c}).maxStatSumPos(1,i) = cluster_shuffled.maxStatSumPos;
        clear cluster_shuffled
    end
    
    if strcmp(ClusterStats, 'SumPos')
        maxStats = ClusterInference.(conditions{c}).maxStatSumPos;
        maxStats = sort(maxStats, 'descend');
        CritVal = maxStats(0.05*size(maxStats,2));
        SigClusters = find(ClusterInference.(conditions{c}).clusters_orig.cluster_statSum > CritVal);
    elseif strcmp(ClusterStats, 'ClusterSize')
        maxStats = ClusterInference.(conditions{c}).maxClusterSizes;
        maxStats = sort(maxStats, 'descend');
        CritVal = maxStats(0.05*size(maxStats,2));
        SigClusters = find(ClusterInference.(conditions{c}).clusters_orig.cluster_size > CritVal);
    end
    
    ClusterInference.(conditions{c}).SigTimePoint = ismember(ClusterInference.(conditions{c}).clusters_orig.cluster_timecourse, SigClusters);
    clear maxStats CritVal SigClusters
end
save(fullfile(OutputDir, 'MEG_RSA_CategoryModel_stats.mat'), 'ClusterInference');