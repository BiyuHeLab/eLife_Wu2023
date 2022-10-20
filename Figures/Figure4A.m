% ENTER a model name
Model = 'TwoState'; % (str) either 'Recognition' or 'TwoState'
if strcmp(Model, 'Recognition') ==1
    bar_color = [0.1216 0.4667 0.7059];
    YAxisLocation = 'right';
elseif strcmp(Model, 'TwoState') ==1
    bar_color = [1 0.4980 0.0549];
     YAxisLocation = 'left';
end

S = load(['../Data/Fusion_' Model 'Model_ClusterInference.mat']);

ROIs = {'V1', 'V2', 'V3', 'placeholder1',...
        'loc_face',  'loc_animal', 'loc_house', 'loc_object', 'placeholder2' ...
        'active_IPS_L', 'active_IPS_R',  'active_MFG_L', 'active_MFG_R', 'active_IFJ_R', ...
        'active_aPCC', 'active_aInsula_L', 'active_aInsula_R', 'active_OFC_R', 'placeholder3'...
        'deactive_AG_L', 'deactive_AG_R', 'deactive_mPFC', 'deactive_PCC',...
        'deactive_SFG_L', 'deactive_SFG_R','deactive_STG_L', 'deactive_STG_R'
        };
    
PlotLabels = {'V1', 'V2', 'V3', '',...
    'face', 'animal', 'house', 'object', ''...
     'L IPS', 'R IPS',  'L MFG', 'R MFG', 'R IFJ',...
     'aPCC', 'L aInsula', 'R aInsula','R OFC', ''...
    'L AG', 'R AG', 'mPFC', 'PCC', ...
    'L SFG', 'R SFG', 'L STG','R STG'};
    
times = -0.5:0.0025:2;
%%
first = zeros(length(ROIs), length(Model));
sig_matrix = nan(length(ROIs), length(times));

for r = 1:length(ROIs)    
    if ~ismember(r, [4, 9, 19])
        sig_matrix(r,:) = S.ClusterInference.(Model).(ROIs{r}).SigTimePoint_SumPos;
        sig_matrix(r,1:200) = 0;
    else
        sig_matrix(r,:) = 1;
    end    
end

figure; hold on
ax = gca;
axis ij
ax.XLim = [-0.2 2];
pbaspect([1 1 1])
box off

ax.YLim = [1-1, length(ROIs)+1];
ax.YTick = 1:length(ROIs);
ax.YTickLabel = PlotLabels;
ax.YAxisLocation = YAxisLocation;
xlabel('Time (s) relative to stimulus onset', 'FontSize', 10)

line([0 0], ax.YLim, 'Color', 'k');
for i = 0.1:0.1:1.9
    line([i i], ax.YLim, 'Color', [0.8 0.8 0.8], 'LineStyle', '--');
end

for r = 1:length(ROIs)
    if ~ismember(r, [4, 9, 19])
        plot(times(sig_matrix(r,:)==1), r, 'Marker', 's',...
             'Markersize', 3, 'MarkerFaceColor', bar_color, 'MarkerEdge', 'none');
    else
        plot(times, r*ones(1, length(times)), 'Color', [0.7 0.7 0.7], 'LineStyle', '--', 'LineWidth', 3);
    end
    
end