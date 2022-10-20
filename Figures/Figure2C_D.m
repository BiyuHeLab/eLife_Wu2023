clear; clc;
load('MEG_avgRDM_100Hz.mat')
times = -0.5:0.01:2;
face = 1:5; house = 6:10; object = 11:15; animal = 16:20; 
%%
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