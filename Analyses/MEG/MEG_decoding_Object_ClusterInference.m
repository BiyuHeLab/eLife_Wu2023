clear; clc;
DataDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/decoding_data/';
addpath('/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/Scripts');

SJs = {'AA', 'AL', 'AR', 'BJB', 'CW', 'DJ', 'EC', 'FSM'...
    'JA', 'JC', 'JP', 'JS', 'LS', 'MC', 'NA', 'NC', 'SL', ...
    'SM', 'TK', 'TL'};
%'AC', 'AW', 'NM', 'SF'
analysis = {'1vs2', '1vs3', '1vs4', '2vs3', '2vs4', '3vs4'};
times = -0.5:0.01:2;
n_perms = 5000;
conditions = {'rec', 'unrec'};
measure = 'balanced_accuracy_minus_chance';
rec = load(fullfile(DataDir, 'Main_ObjectDecodingResults_R.mat'), 'accuracy');
unrec = load(fullfile(DataDir, 'Main_ObjectDecodingResults_U.mat'), 'accuracy');
%%  Averaging across all 6 pair-wise decoding accuracies, for each subject
% respectively.
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
        [pval,~,stats] = signrank(DecodingScore.(conditions{c})(1:20,i_time), 0, 'tail','right');
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
            %pvals(i_perm,i_time) = p; hypothesis(i_perm,i_time) = h;
            zval.(conditions{c})(i_perm,i_time) = stats.zval; %signedrank(i_perm,i_time) = stats.signedrank;
            clear stats
        end
        clear tmp flips FlippedData
    end
end

% Determine significant time point for each permutation sample
% perm_pval.rec = nan(n_perms,length(times));
% perm_pval.unrec = nan(n_perms,length(times));
%
for c = 1:length(conditions)
    for i_perm = 1:n_perms
        for t = 1:length(times)
            pval.(conditions{c})(i_perm,t) = ...
                length(find(perm_zval.(conditions{c})(:,t) > perm_zval.(conditions{c})(i_perm,t))) / n_perms;
        end
    end
end
save(fullfile(DataDir, 'MEG_decoding_Object_shuffled.mat'), 'perm_zval', 'perm_pval')
%% Define significant time points by comparing the original z-stats against
% the empirical null distribution
Orig_pval.rec = nan(1,length(times));
Orig_pval.unrec = nan(1,length(times));
for c = 1:length(conditions)
    for t = 1:length(times)
        Orig_pval.(conditions{c})(1,t) =  ...
            length(find(perm_zval.(conditions{c})(:,t) > Orig_zval.(conditions{c})(1,t))) / n_perms;
    end
end

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
ClusterStats = 'SumPos'; % 'ClusterSize', or 'SumPos'
for c = 1:length(conditions)
    if strcmp(ClusterStats, 'ClusterSize')
        maxStats = ClusterInference.maxClusterSizes.(conditions{c});
        maxStats = sort(maxStats, 'descend');
        CritVal = maxStats(0.05*size(maxStats,2));
        SigClusters.(conditions{c}) = find(ClusterInference.clusters_orig.(conditions{c}).cluster_size > CritVal);
    elseif strcmp(ClusterStats, 'SumPos')
        maxStats = ClusterInference.maxStatSumPos.(conditions{c});
        maxStats = sort(maxStats, 'descend');
        CritVal = maxStats(0.05*size(maxStats,2));
        SigClusters.(conditions{c}) = find(ClusterInference.clusters_orig.(conditions{c}).cluster_statSum > CritVal);
    end
    ClusterInference.maxStatSumPos.SigTimePoint.(conditions{c}) = ...
        ismember(ClusterInference.clusters_orig.(conditions{c}).cluster_timecourse,...
        SigClusters.(conditions{c}));
    clear maxClusterSizes CritVal SigClusters
end
clear Orig_zval Orig_pval perm_zval perm pval
%%
Colors = {
    [0.1725 0.6275 0.1725];% green
    [0.8392 0.1529 0.1569]}; %red
BarPos = [7, 6.8]; 


figure(ceil(100*rand(1)))
% Plot group mean decoding accuracy and SEM
for c = 1:length(conditions)
    mu = nanmean(DecodingScore.(conditions{c}),1);
    sem = nanstd(DecodingScore.(conditions{c}))./sqrt(size(DecodingScore.(conditions{c}),1));
    max_score = max(mu);
    peak.(conditions{c}) = times(find(mu==max_score));
    shadedErrorBar(times,mu, sem, 'lineprops', {'color', Colors{c}, 'LineWidth', 2});
    hold on
    
end

for c = 1:length(conditions)
    if any(ClusterInference.maxStatSumPos.SigTimePoint.(conditions{c})==1)
        sig_time = times(ClusterInference.maxStatSumPos.SigTimePoint.(conditions{c})==1);
        plot(sig_time, BarPos(c), 'Marker', 's',...
            'Markersize', 4, 'MarkerFaceColor', Colors{1}, 'MarkerEdge', 'none')   
    end
end

ax = gca;
ax.XLim = [-0.5 2];
ax.XTick = -0.5:0.5:2;
ax.YLim = [-4.5 7];
ax.YTick = -4:2:7;  
plot(times, zeros(1, 251), 'LineStyle', '--', 'color', 'k');
ax.YTickLabel = {'46', '48', '50', '52', '54', '56'};
line([0 0], ax.YLim, 'color', 'k') 
pbaspect([2 1.3 1])

line([sig_time(1) sig_time(1)],ax.YLim, 'LineStyle', '--', 'color', 'k')
legend('Rcognized', 'Unrecognized') 
% 
xlabel('Time (sec) relative to stimulus onset', 'FontSize', 8, 'Fontweight', 'normal')
ylabel('Decoding Accuracy (%)', 'FontSize', 8, 'Fontweight', 'normal')