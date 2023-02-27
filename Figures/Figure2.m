%PLOT FIGURE 2 A-E

addpath('./helper functions')
load('./Data/MEG_decoding_Recognition.mat', 'accuracy')
load('./Data/MEG_decoding_Recognition_stats.mat')
times = -0.5:0.01:2;
%% 2A Recognition outcome decoding
Colors = [0 0.2 0.4];  
mu = nanmean(accuracy,1);
sem = nanstd(accuracy,1)./sqrt(size(accuracy,1));

figure;hold on
shadedErrorBar(times,smooth(mu,3), smooth(sem,3), 'lineprops', {'color', ...
    Colors, 'LineWidth', 1});

if any(times(ClusterInference.SigTimePoint==1))
    sig_time = times(ClusterInference.SigTimePoint==1);
    plot(sig_time, 9, 'Marker', 's',...
        'Markersize', 4, 'MarkerFaceColor', Colors, 'MarkerEdge', 'none');
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

%% 2B Category decoding
rec = load('./Data/MEG_decoding_Object_R.mat', 'accuracy');
unrec = load('./Data/MEG_decoding_Object_U.mat', 'accuracy');
load('./Data/MEG_decoding_Object_stats.mat', 'ClusterInference');

times = -0.5:0.01:2;
analysis = {'1vs2', '1vs3', '1vs4', '2vs3', '2vs4', '3vs4'};
conditions = {'rec', 'unrec'};
% Averaging across all 6 pair-wise decoding accuracies, for each subject
% respectively.
DecodingScore.rec = nan(size(rec.accuracy.class_1vs2,1), length(times));
DecodingScore.unrec= nan(size(unrec.accuracy.class_1vs2,1), length(times));

for s = 1:size(rec.accuracy.class_1vs2,1)
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

Colors = {[0.1725 0.6275 0.1725],...
    [0.8382, 0.1529, 0.1569]};  
figure;
for i = 1:2
    mu = nanmean(DecodingScore.(conditions{i}),1);
    sem = nanstd(DecodingScore.(conditions{i}),1)./sqrt(size(DecodingScore.(conditions{i}),1));
    
    hold on
    shadedErrorBar(times, mu, sem, 'lineprops', {'color',...
         Colors{i}, 'LineWidth', 1});
    
    if any(times(ClusterInference.maxStatSumPos.SigTimePoint.(conditions{i})==1))
        sig_time = times(ClusterInference.maxStatSumPos.SigTimePoint.(conditions{i})==1);
        plot(sig_time, 7, 'Marker', 's',...
            'Markersize', 4, 'MarkerFaceColor', Colors{i}, 'MarkerEdge', 'none');
    end
end

ax = gca;
ax.XLim = [-0.5 2];
ax.YLim = [-4 7];
ax.YTick = -4:2:15;
plot(times, zeros(1, length(times)), 'LineStyle', '--', 'color', 'k');
line([sig_time(1) sig_time(1)],ax.YLim, 'LineStyle', '--', 'color', 'k')
ax.YTickLabel = {'46','48', '50', '52', '54', '56'};
line([0 0], ax.YLim, 'color', 'k')
pbaspect([2 1.3 1])

xlabel('Time (sec) relative to stimulus onset', 'FontSize', 10, 'Fontweight', 'normal')
ylabel('Decoding accuracy (%)', 'FontSize', 10, 'Fontweight', 'normal')
clear
%% 2C RDMs at different latencies 
load('./Data/MEG_RDM.mat', 'AvgData')
map = viridis;

times = -0.5:0.01:2;
face = 1:5; house = 6:10; object = 11:15; animal = 16:20; 
%
figure(ceil(100*rand(1)));
latency = [0 0.2, 0.29, 0.47] ;
for i = 1:length(latency)
    t = find(abs(times-latency(i)) < 0.0001);
    
    avgRDM = squareform(AvgData(t,:));
    subplot(1,4,i)
    imagesc(avgRDM);
    title([num2str(times(t)*1000) ' ms'], 'FontWeight', 'normal', 'FontSize', 10)
    pbaspect([1 1 1])
    colormap(map)
    ax = gca;
    ax.XTickLabel = '';
    ax.YTickLabel = '';
    ax.CLim = [0 1];
    set(gca, 'XColor', 'none')
    set(gca, 'YColor', 'none')
    grid on;
    ax.XTick = 0.5:5:39.5;
    ax.YTick = 0.5:5:39.5;
    set(gca, 'PlotBoxAspectRatio', [1 1 1])
end
%% 2D Multidimensional Scaling (MDS) at different latencies 
Colors = {[0.1216 0.4687 0.7059];
    [1 0.4980 0.0549];
    [0.1725 0.6275 0.1725];
    [0.8382, 0.1529, 0.1569]};

figure(ceil(100*rand(1)));
for i = 1:length(latency)
        t = find(abs(times-latency(i)) < 0.0001);
      
        avgRDM = squareform(AvgData(t,:));
        Y = mdscale(avgRDM, 2);
      
        subplot(1,4,i)
        %seen stimuli
        scatter(Y(face,1), Y(face,2), 'o', 'MarkerFaceColor', Colors{1}, 'MarkerEdgeColor', 'none'); %blue
        hold on
        scatter(Y(animal,1), Y(animal,2), 'o', 'MarkerFaceColor', Colors{4}, 'MarkerEdgeColor', 'none') %red
        scatter(Y(house,1), Y(house,2), 'o', 'MarkerFaceColor', Colors{2}, 'MarkerEdgeColor', 'none') %orange
        scatter(Y(object,1), Y(object,2), 'o', 'MarkerFaceColor', Colors{3}, 'MarkerEdgeColor', 'none') %green
        %unseen stmuli
        scatter(Y(face+20,1), Y(face+20,2), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Colors{1});
        scatter(Y(animal+20,1), Y(animal+20,2), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Colors{4})
        scatter(Y(house+20,1), Y(house+20,2), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Colors{2})
        scatter(Y(object+20,1), Y(object+20,2), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', Colors{3})
    
    set(gca, 'PlotBoxAspectRatio', [1 1 1])
    set(gca,'XLim',[-0.8,0.8])
    set(gca,'YLim',[-0.8 0.8])
    set(gca,'XTick', [])
    set(gca,'YTick', [])
    box on
    title([num2str(1000*times(t)) ' ms'])

    clear avgRDM Y
end
%clear
%% 2E Across-image dissimilarity for recognized and unrecognized trials
clear; clc;

%addpath('./helper functions')
load('./Data/MEG_RDM.mat', 'Subjects')

times = -0.5:0.01:2;
conditions = {'seen', 'unseen', 'difference'};

SubRDMs.seen = zeros(length(times), size(Subjects,3));
SubRDMs.unseen = zeros(length(times), size(Subjects,3));

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
          [0 0.2 0.4]}; %red-orange
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
%% 2F
