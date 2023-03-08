clear; clc;
% =========================================================================
% Xompute the mean between-exemplar dissimilarity in recognized and
% unrecognized trials, respectively.
% Compute group statistics

% input data:
% Subject's MEG RDV in n_timepoints x n_conditionpairs x n_subjects format
% n_conditionpairs = 780 which corresponds to a 40x40 RDM

% outputs:
% MEG_RSA_MeanDissimilarity.mat
% MEG_RSA_MeanDissimilarity_stats.mat (Stats correspond to the cluster-based
% permutation test)

% Last update: 3/7/2023 @Yuan-hao Wu
% =========================================================================

SubjData = load('ADD PATH TO THE REQUIRED INPUT DATA');
OutputDir = 'PATH FOR THE OUTPUT FILES';
addpath('../../Figures/helper functions')

n_subjects=24;
times = -0.5:0.01:2;
conditions = {'seen', 'unseen', 'difference'};

%%
SubRDMs.seen = zeros(n_subjects, length(times));
SubRDMs.unseen = zeros(n_subjects, length(times));

for sub = 1:n_subjects
    for t_idx = 1:length(times)
        tmp = squareform(SubjData(t_idx,:,sub));
        SubRDMs.seen(sub,t_idx) = nanmean(squareform(tmp(1:20,1:20)));
        SubRDMs.unseen(sub, t_idx) = nanmean(squareform(tmp(21:40,21:40)));
        clear tmp
    end
    clear t_idx
end
save(fullfile(OutputDir, 'MEG_RSA_MeanDissimilarity.mat'), 'SubRDMs') 
%%
for i_time = 1:length(times)
    
    [p,h,stats] = signrank(SubRDMs.seen(:,i_time), 1, 'tail','left');
    Wilcoxon.seen.pval(1,i_time) = p; Wilcoxon.seen.h(1,i_time) = h;
    Wilcoxon.seen.zval(1,i_time) = stats.zval;
    clear p h stats
    
    [p,h,stats] = signrank(SubRDMs.unseen(:,i_time), 1, 'tail','left');
    Wilcoxon.unseen.pval(1,i_time) = p; Wilcoxon.unseen.h(1,i_time) = h;
    Wilcoxon.unseen.zval(1,i_time) = stats.zval;
    clear p h stats
    
    [p,h,stats] = signrank(SubRDMs.seen(:,i_time), SubRDMs.unseen(:,i_time), 'tail','left');
    Wilcoxon.difference.pval(1,i_time) = p; Wilcoxon.difference.h(1,i_time) = h;
    Wilcoxon.difference.zval(1,i_time) = stats.zval; 
    clear p h stats
end
%% SIGN PERMUTATION TEST
% 1.) randomly flips the signs of subjects' data with 50%
% 2.) run Wilcoxon test with the data
% 3.) repeat step 1 and 2 for 5000 times
rng('default')
rng(12332)
n_perms = 5000;

for i_perm = 1:n_perms
    % sample a a random subset of subjects whose data will be sign flipped
    n_flips = randsample(1:n_subjects,1);  
    flipped_samples = randsample(1:n_subjects, n_flips, false);
    
    % now flip the sign of these subjects' data
    % Pearson's distance has an expected value of 1. The flipped value of
    % a distance of 1.1 would be 0.9, not -1.1.
    % Compute the flipped value using 1 + 1-Pearson's distance)
    FlippedData = SubRDMs;
    
    for c=1:2
        FlippedData.(conditions{c})(flipped_samples,:) = 1+ (1-FlippedData.(conditions{c})(flipped_samples,:)); %.*-1;
    end
   
    for i_time = 1:length(times)
        [pval,h,stats] = signrank(FlippedData.(conditions{1})(:,i_time), 1, 'tail','left');
        PermStats.(conditions{1}).pval(i_perm,i_time) = pval;
        PermStats.(conditions{1}).h(i_perm,i_time) = h;
        PermStats.(conditions{1}).zval(i_perm,i_time)  = stats.zval;
        clear pval stats h
        
        [pval,h,stats] = signrank(FlippedData.(conditions{2})(:,i_time), 1, 'tail','left');
        PermStats.(conditions{2}).pval(i_perm,i_time) = pval;
        PermStats.(conditions{2}).h(i_perm,i_time) = h;
        PermStats.(conditions{2}).zval(i_perm,i_time)  = stats.zval;
        clear pval stats h
        
        [pval,h,stats] = signrank(FlippedData.(conditions{1})(:,i_time),...
            FlippedData.(conditions{2})(:,i_time), 'tail','left');
        PermStats.(conditions{3}).pval(i_perm,i_time) = pval;
        PermStats.(conditions{3}).h(i_perm,i_time) = h;
        PermStats.(conditions{3}).zval(i_perm,i_time)  = stats.zval;
        clear pval stats h   
    end 
    clear flipped_samples n_flips FlippedData
end
%% cluster inference
CDT = 0.05;  
ClusterStats = 'SumAbs';  % ['SumPos', 'ClusterSize]
conditions = {'seen', 'unseen', 'difference'};

for c = 1:length(conditions)
    clusters_orig = find_temporal_clusters(Wilcoxon.(conditions{c}).zval , Wilcoxon.(conditions{c}).pval, CDT);
    ClusterInference.(conditions{c}).clusters_orig = clusters_orig;
    clear clusters_orig
    
    for i = 1:n_perms
        cluster_shuffled = find_temporal_clusters(PermStats.(conditions{c}).zval(i,:),PermStats.(conditions{c}).pval(i,:), CDT);
        
        ClusterInference.(conditions{c}).maxClusterSizes(1,i) = cluster_shuffled.maxSize;
        ClusterInference.(conditions{c}).maxStatSumAbs(1,i) = cluster_shuffled.maxStatSumAbs;
        clear cluster_shuffled
    end
    
    if strcmp(ClusterStats, 'SumAbs')
        maxStats = ClusterInference.(conditions{c}).maxStatSumAbs;
        maxStats = sort(maxStats, 'descend');
        CritVal = maxStats(0.05*size(maxStats,2));
        SigClusters = find(abs(ClusterInference.(conditions{c}).clusters_orig.cluster_statSum) > CritVal);
    elseif strcmp(ClusterStats, 'ClusterSize')
        maxStats = ClusterInference.(conditions{c}).maxClusterSizes;
        maxStats = sort(maxStats, 'descend');
        CritVal = maxStats(0.05*size(maxStats,2));
        SigClusters = find(ClusterInference.(conditions{c}).clusters_orig.cluster_size > CritVal);
    end
    
    ClusterInference.(conditions{c}).SigTimePoint = ismember(ClusterInference.(conditions{c}).clusters_orig.cluster_timecourse, SigClusters);
    clear maxStats CritVal SigClusters
end
save(fullfile(OutputDir, 'MEG_RSA_MeanDissimilarity_stats.mat'), 'ClusterInference', 'PermStats') 