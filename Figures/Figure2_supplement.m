clear; clc
load('../Data/MEG_decoding_HighContrast.mat')
load('../Data/MEG_decoding_HighContrast_shuffled.mat')
analysis = {'1vs2', '1vs3', '1vs4', '2vs3', '2vs4', '3vs4'};
times = -0.1:0.01:0.8;
n_perms = 5000; 

%% compute mean accuracy across all pairwise comparisons
% and Wicoxon sign rank test
mean_accuracy = nan(size(accuracy.class_1vs2,1), length(times));
for s = 1:size(accuracy.class_1vs2,1)
    tmp1 = [];
    for i_pair = 1:length(analysis)
        tmp1 = [tmp1; accuracy.(['class_' analysis{i_pair}])(s,:)];
    end
    mean_accuracy(s,:) = nanmean(tmp1,1);
end
clear tmp1

for i_time = 1:length(times)
    [p,~,stats] = signrank(mean_accuracy(:,i_time), 0, 'tail','right');
    Orig_pval(1,i_time) = p; 
    Orig_zval(1,i_time) = stats.zval;
    clear p h stats
end
clear accuracy analysis i_pair measures

%% Define significant time points by comparing the original z-stats against

% find clusters of significant time points in the original sample
clusters_orig = find_temporal_clusters(Orig_zval(1,:),Orig_pval(1,:), 0.05);
ClusterInference.clusters_orig = clusters_orig;
clear clusters_orig

% generate maxClusterSize and max SumStat null distribution
for i_perm = 1:n_perms
    cluster_shuffled = find_temporal_clusters(perm_zval(i_perm,:), perm_pval(i_perm,:), 0.05);
    ClusterInference.maxClusterSizes(1,i_perm) = cluster_shuffled.maxSize;
    ClusterInference.masStatSumPos(1,i_perm) = cluster_shuffled.maxStatSumPos;
    clear cluster_shuffled
end
%%
maxStats = ClusterInference.maxClusterSizes;
maxStats = sort(maxStats, 'descend');
CritVal = maxStats(0.05*size(maxStats,2));
SigClusters = find(ClusterInference.clusters_orig.cluster_size > CritVal);

ClusterInference.SigTimePoint = ismember(ClusterInference.clusters_orig.cluster_timecourse, SigClusters);
clear maxClusterSizes CritVal SigClusters

%% PLOT decoding accuracies over time

Colors = {
    [0.1725 0.6275 0.1725];% green
 }; %red

figure(ceil(100*rand(1))) 
mu = nanmean(mean_accuracy,1);
sem = nanstd(mean_accuracy)./sqrt(size(mean_accuracy,1));
y2 = shadedErrorBar(times,smooth(mu(1,:),1), sem(1,:), 'lineprops',...
    {'color', [0.1725 0.6275 0.1725], 'LineWidth', 2});
hold on

if ~isempty(times(ClusterInference.SigTimePoint==1))
    sig_time = times(ClusterInference.SigTimePoint==1);
    plot(sig_time, 15, 'Marker', 's', 'Markersize', 4, 'MarkerFaceColor',...
        [0.1725 0.6275 0.1725], 'MarkerEdge', 'none');
end

ax = gca;
ax.XLim = [-0.1 0.8];
ax.XTick = -0.1:0.1:0.8;
ax.YLim = [-2 15];
ax.YTick = [-2, 0, 5 10 15];  
plot(times, zeros(1, 91), 'LineStyle', '--', 'color', 'k');
ax.YTickLabel = {'48', '50', '55', '60', '65'};
line([0 0], ax.YLim, 'color', 'k')
line([sig_time(1) sig_time(1)],ax.YLim, 'LineStyle', '--', 'color', 'k')
xlabel('Time (sec) relative to stimulus onset', 'FontSize', 8, 'Fontweight', 'normal')
ylabel('Decoding Accuracy (%)', 'FontSize', 8, 'Fontweight', 'normal')
pbaspect([2 1.3 1])
clear 