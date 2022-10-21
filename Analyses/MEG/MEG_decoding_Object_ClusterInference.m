clear; clc;
DataDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/decoding_data/';
addpath('/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/Scripts');

SJs = {'AA', 'AC','AL', 'AR', 'AW', 'BJB', 'CW', 'DJ', 'EC', 'FSM'...
    'JA', 'JC', 'JP', 'JS', 'LS', 'MC', 'NA', 'NC', 'NM', 'SF ', 'SL', ...
    'SM', 'TK', 'TL'};
analysis = {'1vs2', '1vs3', '1vs4', '2vs3', '2vs4', '3vs4'};
measure = 'balanced_accuracy_minus_chance';
times = -0.5:0.01:2;
n_perms = 5000;
conditions = {'rec', 'unrec'};
%% concatenate decoding accuracies yielded from all time points
for i= 1:length(analysis)
    for s = 1:length(SJs)
        SubDir = fullfile(DataDir, 'RecObj', SJs{s}, ['s_ensemble_trial_5folds_' analysis{i}]);
        timecourses = nan(1,length(times));
        for i_time = 1:length(times)
            curr_Dir = fullfile(SubDir, ['t' num2str(i_time, '%03.f')]);
            load(fullfile(curr_Dir, ['res_' measure '.mat']));
            timecourses(1,i_time) = results.(measure).output;
            clear results curr_Dir
        end
        accuracy.(['class_' analysis{i}])(s,:) = timecourses(1,:);
        clear timecourses SubDir i_time k
    end
end
save(fullfile(DataDir, 'MEG_decoding_Object_R.mat'), 'accuracy')

for i= 1:length(analysis)
    for s = 1:length(SJs)
        SubDir = fullfile(DataDir, 'UnrecObj', SJs{s}, ['s_ensemble_trial_5folds_' analysis{i}]);
        timecourses = nan(1,length(times));
        for i_time = 1:length(times)
            curr_Dir = fullfile(SubDir, ['t' num2str(i_time, '%03.f')]);
            load(fullfile(curr_Dir, ['res_' measure '.mat']));
            timecourses(1,i_time) = results.(measure).output;
            clear results curr_Dir
        end
        accuracy.(['class_' analysis{i}])(s,:) = timecourses(1,:);
        clear timecourses SubDir i_time k
    end
end
save(fullfile(DataDir, 'MEG_decoding_Object_U.mat'), 'accuracy')
%%  Averaging across all 6 pair-wise decoding accuracies, for each subject
% respectively.
rec = load(fullfile(DataDir, 'MEG_decoding_Object_R.mat'), 'accuracy');
unrec = load(fullfile(DataDir, 'MEG_decoding_Object_U.mat'), 'accuracy');

DecodingScore.rec = nan(length(SJs), length(times));
DecodingScore.unrec= nan(length(SJs), length(times));


for s = 1:length(SJs)
    tmp1 = [];
    tmp2 = [];
    for i_pair = 1:length(analysis)
        tmp1 = [tmp1; rec.accuracy.(['class_' analysis{i_pair}])(s,:)];
        tmp2 = [tmp2; unrec.accuracy.(['class_' analysis{i_pair}])(s,:)];
    end
    DecodingScore.rec(s,:) = nanmean(tmp1,1);
    DecodingScore.unrec(s,:) = nanmean(tmp2,1);
end
clear tmp1 tmp2 rec unrec
%% Wicoxon sign rank test on observed data
for c = 1:length(conditions)
    for i_time = 1:length(times)
        [pval,~,stats] = signrank(DecodingScore.(conditions{c})(:,i_time), 0, 'tail','right');
        Orig_pval.(conditions{c})(1,i_time) = pval;
        Orig_zval.(conditions{c})(1,i_time) = stats.zval;
        clear p h stats
    end
end
clear i_pair i_time analysis
%% Run flip-sign test and generate max stats distribution
% 1.) flip the sign of a random subset of subject data
% 2.) run statistical test with the data
% 3.) repeat step 1 and 2 for 500 times
% 4.) generate a max sumstat distribution
rng('default')
rng(12332)
for c = 1:length(conditions)
    for i_perm = 1:n_perms
        tmp = randsample(1:size(SJs,2),1);  % number of subjects whose data will be sign flipped
        flips = randsample(1:size(SJs,2), tmp); % randomly choose n subjects
        % now flip the sign of these subjects' data
        FlippedData = DecodingScore.(conditions{c});
        FlippedData(flips,:) = FlippedData(flips,:).*-1;
        % run Wilcoxon test on the the dataset including sign-flipped subjects
        for i_time = 1:length(times)
            [~,~,stats] = signrank(FlippedData(:,i_time), 0, 'tail','right');
            pvals(i_perm,i_time) = p;
            zval.(conditions{c})(i_perm,i_time) = stats.zval; %signedrank(i_perm,i_time) = stats.signedrank;
            clear stats
        end
        clear tmp flips FlippedData
    end
end

save(fullfile(DataDir, 'MEG_decoding_Object_shuffled.mat'), 'perm_zval', 'perm_pval')
%% Define significant time points by comparing the original z-stats against
% the empirical null distribution

% find clusters of significant time points in the original sample
for c = 1:length(conditions)
    clusters_orig.(conditions{c}) = ...
        find_temporal_clusters(Orig_zval.(conditions{c})(1,:),Orig_pval.(conditions{c})(1,:), 0.05);
end
ClusterInference.clusters_orig = clusters_orig;
clear clusters_orig

% generate maxClusterSize and max SumStat null distribution
for c = 1:length(conditions)
    for i_perm = 1:n_perms
        cluster_shuffled = find_temporal_clusters(perm_zval.(conditions{c})(i_perm,:),...
            perm_pval.(conditions{c})(i_perm,:), 0.05);
        ClusterInference.maxClusterSizes.(conditions{c})(1,i_perm) = cluster_shuffled.maxSize;
        ClusterInference.maxStatSumPos.(conditions{c})(1,i_perm) = cluster_shuffled.maxStatSumPos;
        clear cluster_shuffled
    end
end

%%
for c = 1:length(conditions)
    maxStats = ClusterInference.maxStatSumPos.(conditions{c});
    maxStats = sort(maxStats, 'descend');
    CritVal = maxStats(0.05*size(maxStats,2));
    SigClusters.(conditions{c}) = find(ClusterInference.clusters_orig.(conditions{c}).cluster_statSum > CritVal);
    
    ClusterInference.maxStatSumPos.SigTimePoint.(conditions{c}) = ...
        ismember(ClusterInference.clusters_orig.(conditions{c}).cluster_timecourse,...
        SigClusters.(conditions{c}));
    clear maxClusterSizes CritVal SigClusters
end
save(fullfile(pwd, 'MEG_decoding_Object_stats.mat'), 'ClusterInference')
clear