% =========================================================================
% created on: unkown
% Plots the commonality coefficient of MEG,, fMRI  and a single model with
% at each time point between -0.5 and 2.0, for each ROI separately.
% Correspond to Figure 4B in the Manuscript
% The commonality coefficient reflects the part of the MRI-MEG covariance
% that can be attributed to a given model.
%author: @Yuan-hao Wu
% last update: 2022/09/27 @Yuan-hao Wu
% =========================================================================

clear; clc
addpath('./helper functions')
Models = {'Recognition', 'TwoState'};

load('./Data/Fusion_Commonality.mat');
load('./Data/Fusion_MEG-MRI_Correlation.mat');

ROIs = {'V1', 'V2', 'V3', 'loc_face', 'loc_animal', 'loc_house', 'loc_object',  ...
    'active_IPS_L', ...
    'active_IPS_R', 'active_aPCC', 'active_aInsula_L', 'active_aInsula_R', ...
    'active_IFJ_R', 'active_MFG_L', 'active_MFG_R', 'active_OFC_R',...
    'deactive_AG_L', 'deactive_AG_R', ...
    'deactive_mPFC', 'deactive_PCC', 'deactive_SFG_L', 'deactive_SFG_R', ...
    'deactive_STG_L', 'deactive_STG_R'};

PlotLabels = {'V1', 'V2', 'V3', 'face', 'animal', 'house', 'object'...
    'L IPS', 'R IPS', 'aPCC', 'L aInsula', 'R aInsula', 'R IFJ',...
    'L MFG', 'R MFG', 'R OFC',...
    'L AG', 'R AG', 'mPFC', 'PCC', 'L SFG', 'R SFG', 'L STG', 'R STG'};
%% Plot with quadratic scale
times = -0.5:0.0025:2;
colors = {[0.1216 0.4667 0.7059], [1 0.4980 0.0549]}; 
tytick = [-0.08 -0.04 -0.02:0.01:0.02 0.04 0.16 0.20];
tyticklabel = num2cell(tytick);
tylim = tytick([1 end]);
bar_pos = [0.19, 0.17];

 for r = 1:length(PlotLabels)
    h = figure(r);
    hold on
    
    ha = gca;
    ha.YTickMode = 'manual';
    ha.YTick = sign(tytick).*sqrt(abs(tytick));
    ha.YTickLabel = tyticklabel;
    ha.YLim = sign(tylim).*sqrt(abs(tylim));
    ha.XLim = [min(times) max(times)];
    pbaspect([2 1.5 1])
   
    ylabel('Commonality coefficient', 'FontSize', 10, 'Fontweight', 'normal')
    xlabel('Time (s) relative to stimulus onset', 'FontSize', 10, 'Fontweight', 'normal')
    title(PlotLabels{r})
        
    % plot MEG-fMRI fusion as shaded gray area
    UpperB = AvgRho.(ROIs{r}).^2;
    h1 = area(times, sqrt(UpperB));
    h1.FaceColor = [0.7843 0.7843 0.7843];
    h1.LineStyle = 'none';
    
    % plot commonality coefficient for each model
    for m = 1:length(Models)        
        Rsquared_orig = C.(Models{m}).(ROIs{r});
        s = sign(Rsquared_orig);
        h2 = plot(times, s.*sqrt(abs(Rsquared_orig)), 'LineWidth', 2, 'Color', colors{m});
        clear Rsquared_orig h2 s
      
    % add horizontal significance bars and onset times
        load(['./Data/Fusion_' Models{m} 'Model_ClusterInference.mat']);      
        
        if ~isempty(times(ClusterInference.(Models{m}).(ROIs{r}).SigTimePoint_SumPos==1))
            masked_time = times(ClusterInference.(Models{m}).(ROIs{r}).SigTimePoint_SumPos==1);          
            plot(masked_time(masked_time>0), sign(bar_pos(m))*sqrt(abs(bar_pos(m))), 'Marker', 's',...
                'Markersize', 3, 'MarkerFaceColor', colors{m}, 'MarkerEdge', 'none');
            onsets =  masked_time(find(masked_time > 0, 1));
            if ~isempty(onsets)
                line([onsets onsets], ha.YLim, 'LineWidth', 0.5, 'LineStyle', '--', 'Color', colors{m})               
            end
        end        
    end
    line([0 0], ha.YLim, 'Color', 'k');
    plot(times, zeros(1, length(times)), 'Color', 'k', 'LineStyle', '--');   
    clear h h1 ha UpperB 
end
clear