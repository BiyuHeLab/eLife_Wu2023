clear; clc;

addpath('./helper functions')
load('./Data/MEG_RDM.mat', 'Subjects')

times = -0.5:0.01:2;
conditions = {'seen', 'unseen', 'difference'};

SubRDMs.seen = zeros(length(times), size(Subjects,3));
SubRDMs.unseen = zeros(length(times), size(Subjects,3));
%%
for sub = 1:size(Subjects,3)
    for t_idx = 1:length(times)
        tmp = squareform(Subjects(t_idx,:,sub));
        SubRDMs.seen(t_idx, sub) = nanmean(squareform(tmp(1:20,1:20)));
        SubRDMs.unseen(t_idx, sub) = nanmean(squareform(tmp(21:40,21:40)));
        SubRDMs.difference(t_idx, sub) = SubRDMs.seen(t_idx, sub) - SubRDMs.unseen(t_idx, sub);
        clear tmp
    end
    clear t_idx
end
%%
for i_time = 1:length(times)
    
    [p,h,stats] = signrank(SubRDMs.seen(i_time,:), 1, 'tail','left');
    Wilcoxon.seen.p(1,i_time) = p; Wilcoxon.seen.h(1,i_time) = h;
    Wilcoxon.seen.zval(1,i_time) = stats.zval; Wilcoxon.seen.signedrank(1,i_time) = stats.signedrank;
    clear p h stats
    
    [p,h,stats] = signrank(SubRDMs.unseen(i_time,:), 1, 'tail','left');
    Wilcoxon.unseen.p(1,i_time) = p; Wilcoxon.unseen.h(1,i_time) = h;
    Wilcoxon.unseen.zval(1,i_time) = stats.zval; Wilcoxon.unseen.signedrank(1,i_time) = stats.signedrank;
    clear p h stats
    
    [p,h,stats] = signrank(SubRDMs.difference(i_time,:), 0, 'tail','left');
    Wilcoxon.difference.p(1,i_time) = p; Wilcoxon.difference.h(1,i_time) = h;
    Wilcoxon.difference.zval(1,i_time) = stats.zval; Wilcoxon.difference.signedrank(1,i_time) = stats.signedrank;
    clear p h stats
end
%%
colors = {[0.1725 0.6275 0.1725]; %green
    [0.8382, 0.1529, 0.1569];% red
    [0.8500 0.3250 0.0980]}; %red-orange
BarPos = [0.5, 0.47, 0.44];

figure(round(100*rand(1)))
for c = 1:length(conditions)
    mu = nanmean(SubRDMs.(conditions{c}),2);
    SE = nanstd(SubRDMs.(conditions{c})')./sqrt(size(Subjects,3));
    if c==1
        y1 = shadedErrorBar(times, mu, SE, 'lineprops', {'color', colors{c}, 'LineWidth', 2});      
    elseif c==2
         y2 = shadedErrorBar(times, mu, SE, 'lineprops', {'color', colors{c}, 'LineWidth', 2});   
    end
     hold on
%     if ~isempty(times(ClusterInference.(conditions{c}).SigTimePoint))
%         plot(times(ClusterInference.(conditions{c}).SigTimePoint==1), BarPos(c), 'Marker', 's',...
%         'Markersize', 4, 'MarkerFaceColor', colors{c}, 'MarkerEdge', 'none');
%     end


     if ~isempty(times(Wilcoxon.(conditions{c}).h==1))
         plot(times(Wilcoxon.(conditions{c}).h==1), BarPos(c), 'Marker', 's',...
             'Markersize', 4, 'MarkerFaceColor', colors{c}, 'MarkerEdge', 'none');
     end    
end
    
pbaspect([2 1.3 1])
ax= gca;
ax.YLim = [0.4, 1.1];
ax.YTick = ax.YLim(1):0.2:ax.YLim(2);
plot(times, ones(1, length(times)), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
line([0 0],ax.YLim,  'LineWidth', 1, 'Color', 'k')
xlabel('time (sec)')
ylabel('1 - Pearson''s r')