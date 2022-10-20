%
load('MEG_decoding_Recognition.mat', 'accuracy')
load('MEG_decoding_Recognition_shuffled.mat', 'perm_signedrank', 'perm_pval1')

times = -0.5:0.01:2;
n_perms = 5000;

%% Run Wilcoxon test on original data
accuracy=accuracy(1:size(accuracy,1),:);
for i_time = 1:length(times)
    [pval,~,stats] = signrank(accuracy(:,i_time), 0, 'tail','right');
     Orig_signedrank(1,i_time) = stats.signedrank;
     Orig_pval(1,i_time) = pval;
     clear p stats
end

clusters_orig = find_temporal_clusters(Orig_signedrank(1,:),Orig_pval(1,:), 0.05);
ClusterInference.clusters_orig = clusters_orig;
clear clusters_orig

for i_perm = 1:n_perms
    cluster_shuffled = find_temporal_clusters(perm_signedrank(i_perm,:),perm_pval1(i_perm,:), 0.05);
    ClusterInference.maxClusterSizes(1,i_perm) = cluster_shuffled.maxSize;
    ClusterInference.masStatSumPos(1,i_perm) = cluster_shuffled.maxStatSumPos;
    clear cluster_shuffled
end

maxStats = ClusterInference.maxClusterSizes;
maxStats = sort(maxStats, 'descend');
CritVal = maxStats(0.05*size(maxStats,2));
SigClusters = find(ClusterInference.clusters_orig.cluster_size > CritVal);
ClusterInference.SigTimePoint = ismember(ClusterInference.clusters_orig.cluster_timecourse, SigClusters);
clear maxClusterSizes CritVal SigClusters
%%
Colors = [0 0.4470 0.7410];  
mu = nanmean(accuracy,1);
sem = nanstd(accuracy,1)./sqrt(size(accuracy,1));

figure;hold on
shadedErrorBar(times,smooth(mu,3), smooth(sem,3), 'lineprops', {'color',...
    [0 0.4470 0.7410], 'LineWidth', 1});

if any(times(ClusterInference.SigTimePoint==1))
    sig_time = times(ClusterInference.SigTimePoint==1);
    plot(sig_time, 9, 'Marker', 's',...
        'Markersize', 4, 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdge', 'none');
end

ax = gca;
ax.XLim = [-0.5 2];
ax.YLim = [-2 9];
ax.YTick = -2:2:15;
plot(times, zeros(1, length(times)), 'LineStyle', '--', 'color', 'k');
line([sig_time(1) sig_time(1)],ax.YLim, 'LineStyle', '--', 'color', 'k')
ax.YTickLabel = {'', '50', '52', '54', '56', '58'};
line([0 0], ax.YLim, 'color', 'k')
pbaspect([2 1.3 1])

xlabel('Time (sec) relative to stimulus onset', 'FontSize', 10, 'Fontweight', 'normal')
ylabel('Decoding accuracy (%)', 'FontSize', 10, 'Fontweight', 'normal')
clear  