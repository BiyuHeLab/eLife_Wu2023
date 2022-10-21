%PLOT FIGURE 2 A-D


load('./Data/MEG_decoding_Recognition.mat', 'accuracy')
load('./Data/MEG_decoding_Recognition_stats.mat')
times = -0.5:0.01:2;
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

%%
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
%clear


%%
load('./Data/MEG_avgRDM_100Hz.mat')
times = -0.5:0.01:2;
face = 1:5; house = 6:10; object = 11:15; animal = 16:20; 
%
h = figure(ceil(100*rand(1)));
latency = [0 0.2, 0.29, 0.47] ;
for i = 1:length(latency)
    t = find(abs(times-latency(i)) < 0.0001);
    
    avgRDM = squareform(AvgData(t,:));
    subplot(1,4,i)
    imagesc(avgRDM);
    title([num2str(times(t)*1000) ' ms'], 'FontWeight', 'normal', 'FontSize', 10)
    pbaspect([1 1 1])
    %colorbar
    ax = gca;
    ax.XTickLabel = '';
    ax.YTickLabel = '';
    ax.CLim = [0 1];
    set(gca, 'XColor', 'none')
    set(gca, 'YColor', 'none')
    grid on;
    ax.XTick = [0.5:5:39.5];
    ax.YTick = [0.5:5:39.5];
    set(gca, 'PlotBoxAspectRatio', [1 1 1])
end
%%
h = figure(ceil(100*rand(1)));
for i = 1:length(latency)
        t = find(abs(times-latency(i)) < 0.0001);
      
        avgRDM = squareform(AvgData(t,:));
        Y = mdscale(avgRDM, 2);
      
        subplot(1,4,i)
        %seen stimuli
        h = scatter(Y(face,1), Y(face,2), 'o', 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'none'); %red
        hold on
        scatter(Y(animal,1), Y(animal,2), 'o', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', 'none') %orange
        scatter(Y(house,1), Y(house,2), 'o', 'MarkerFaceColor', [0 0.8 0.8], 'MarkerEdgeColor', 'none') %turquioise
        scatter(Y(object,1), Y(object,2), 'o', 'MarkerFaceColor', [0 0 0.8], 'MarkerEdgeColor', 'none') %blue
        %unseen stmuli
        scatter(Y(face+20,1), Y(face+20,2), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [1 0 0]);
        scatter(Y(animal+20,1), Y(animal+20,2), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [1 0.5 0])
        scatter(Y(house+20,1), Y(house+20,2), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0 0.8 0.8])
        scatter(Y(object+20,1), Y(object+20,2), 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0 0 0.8])
    
    set(gca, 'PlotBoxAspectRatio', [1 1 1])
    set(gca,'XLim',[-0.8,0.8])
    set(gca,'YLim',[-0.8 0.8])
    set(gca,'XTick', [])
    set(gca,'YTick', [])
    box on
    title([num2str(1000*times(t)) ' ms'])

    clear avgRDM Y
end
