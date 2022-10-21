%
load('../Data/MEG_decoding_Recognition.mat', 'accuracy')
%load('../Data/MEG_decoding_Recognition_shuffled.mat', 'perm_signedrank', 'perm_pval1')
load('../Data/MEG_decoding_Recognition_ClusterInference.mat')
times = -0.5:0.01:2;
n_perms = 5000;
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