% FIGURE 3 - supplement1-3 
%computes the correlations between the 40x40 group-average fMRI RDM and
% the 40x40 group-average MEG RDM.

clear; clc
load('./Data/Fusion_MEG-MRI_Correlation.mat');
load('./Data/Fusion_MEG-MRI_Correlation_ClusterInference.mat')
ROIs = fieldnames(AvgRho);
PlotLabels = {'V1', 'V2', 'V3', 'V4', 'face', 'animal', 'house', 'object', ...
     'L IPS', 'R IPS',  'aPCC', 'L aInsula', 'R aInsula',  'R IFJ', 'L MFG', 'R MFG','R OFC',...
    'L AG', 'R AG', 'mPFC', 'PCC', 'L SFG', 'R SFG', 'L STG','R STG'};
times = -0.5:0.0025:2;

%%
for r = 1:length(PlotLabels)
    
    figure(r)  
    plot(times, AvgRho.(ROIs{r}), 'LineWidth', 1, 'Color', [0 0.4470 0.7410])
    hold on
  
    if any(times(ClusterInference.(ROIs{r}).SigTimePoint_SumPos==1))
        sig_time = times(ClusterInference.(ROIs{r}).SigTimePoint_SumPos==1);
        plot(sig_time, 0.45, 'Marker', 's',...
            'Markersize', 2, 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerEdge', 'none');
    end
     
    ax = gca;
    pbaspect([1.5 1 1])
    ax.YLim = [-0.2 0.5];
    plot(times, zeros(1, length(times)), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    line([0 0], ylim, 'Color', 'k')
    box off
    ylabel('Spearman''s rho', 'FontSize', 10, 'Fontweight', 'normal')
    xlabel('Time (s) relative to stimulus onset', 'FontSize', 10, 'Fontweight', 'normal')
    title(PlotLabels{r})
end