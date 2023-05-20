clear; clc;
SJs = {'AA', 'AC', 'AL', 'AR', 'AW','BJB', 'CW', 'DJ', 'EC', 'FSM'...
    'JA', 'JC', 'JP', 'JS', 'LS', 'MC', 'NA', 'NC', 'NM', 'SF', 'SL', ...
    'SM', 'TK', 'TL'}; 
%
DataDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/decoding_data/';
measure = 'balanced_accuracy_minus_chance';
times = -0.5:0.01:2;
n_perms = 5000;
%% 
for s = 1:length(SJs)
    
    SubDir = fullfile(DataDir, 'RecognitionState', SJs{s}, 'ensemble_trial_5folds');
    timecourses = nan(1,length(times));
    for i_time = 1:length(times)
        curr_Dir = fullfile(SubDir, ['t' num2str(i_time, '%03.f')]); 
        load(fullfile(curr_Dir, ['res_' measure '.mat']));
        timecourses(1,i_time) = results.(measure).output;
        clear results curr_Dir
    end
    accuracy(s,:) = timecourses(1,:);
    clear timecourses SubDir i_time
    
end
save(fullfile(pwd, 'MEG_decoding_Recognition.mat'), 'accuracy')
%% Run Wilcoxon test on original data
% if determine the CDT by wilcoxon signed rank test
for i_time = 1:length(times)
    [pval,~,stats] = signrank(accuracy(:,i_time), 0, 'tail','right');
     Orig_signedrank(1,i_time) = stats.signedrank;
     Orig_pval(1,i_time) = pval;
     clear p stats
end

% only if determine the CDT by pval yielded from comparing the null distribution 
% Orig_pval = nan(1,length(times));
% for t = 1:length(times)
%     Orig_pval(1,t) = length(find(perm_signedrank(:,t) > Orig_signedrank(1,t))) / n_perms;
% end

% find clusters of significant time points in the original sample
clusters_orig = find_temporal_clusters(Orig_signedrank(1,:),Orig_pval(1,:), 0.05);
ClusterInference.clusters_orig = clusters_orig;
clear clusters_orig
%% SIGN PERMUTATION TEST
% 1.) randomly flips the signs of subjects' data with 50%
% 2.) run Wilcoxon test with the data
% 3.) repeat step 1 and 2 for 5000 times
rng('default')
rng(12332)
perm_pval1= nan(n_perms, length(times));
for i_perm = 1:n_perms
    % sample a a random subset of subjects whose data will be sign flipped
    n_flips = randsample(1:size(accuracy,1),1);  
    flipped_samples = randsample(1:size(accuracy,1), n_flips, false);
    % now flip the sign of these subjects' data
    FlippedData = accuracy; 
    FlippedData(flipped_samples,:) = FlippedData(flipped_samples,:).*-1;
    % run Wilcoxon test on the the dataset including sign-flipped subjects
    for i_time = 1:length(times)
        [pval,~,stats] = signrank(FlippedData(:,i_time), 0, 'tail','right');
        perm_pval1(i_perm,i_time) = pval;
        perm_signedrank(i_perm,i_time)  = stats.signedrank;
        clear pval stats
    end
    clear tmp flips FlippedData
end

% if coompare against null distribution 
perm_pval2= nan(n_perms, length(times));
for i_perm = 1:n_perms
    for t = 1:length(times)
        perm_pval2(i_perm,t) = length(find(perm_signedrank(:,t)...
            > perm_signedrank(i_perm,t))) / n_perms;
    end
end
load(fullfile(pwd, 'MEG_decoding_Recognition_shuffled.mat'), 'perm_pval1', 'perm_pval2', 'perm_signedrank')
%% generate maxClusterSize and max StatSumPos null distribution
for i_perm = 1:n_perms
    cluster_shuffled = find_temporal_clusters(perm_signedrank(i_perm,:),perm_pval1(i_perm,:), 0.05);
    ClusterInference.maxClusterSizes(1,i_perm) = cluster_shuffled.maxSize;
    ClusterInference.masStatSumPos(1,i_perm) = cluster_shuffled.maxStatSumPos;
    clear cluster_shuffled
end

ClusterStats = 'SumPos'; % 'ClusterSize', or 'SumPos'

if strcmp(ClusterStats, 'ClusterSize')
    maxStats = ClusterInference.maxClusterSizes;
    maxStats = sort(maxStats, 'descend');
    CritVal = maxStats(0.05*size(maxStats,2));
    SigClusters = find(ClusterInference.clusters_orig.cluster_size > CritVal);
elseif strcmp(ClusterStats, 'SumPos')
    maxStats = ClusterInference.masStatSumPos;
    maxStats = sort(maxStats, 'descend');
    CritVal = maxStats(0.05*size(maxStats,2));
    SigClusters = find(ClusterInference.clusters_orig.cluster_statSum > CritVal);
end
ClusterInference.SigTimePoint = ismember(ClusterInference.clusters_orig.cluster_timecourse, SigClusters);
clear maxClusterSizes CritVal SigClusters
save(fullfile(DataDir, 'MEG_decoding_Recognition_ClusterInference.mat'))