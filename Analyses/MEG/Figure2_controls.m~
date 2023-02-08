
clear; clc;
addpath('./helper functions')
load('./Data/MEG_RDM.mat', 'Subjects')
map = viridis;
times=-0.5:0.01:2;
n = size(Subjects,3);
n_conditions = 40;

%% Plot individual subjects' MEG RDM at a given time point 
latency = 0.47; %[0 0.2, 0.29, 0.47] ;
t = find(abs(times-latency(1)) < 0.0001);

figure
for sub=1:n
  
    subplot(4,6,sub)
    sub_rdm = squareform(Subjects(t,:,sub));
    imagesc(sub_rdm);
    pbaspect([1 1 1])

    colormap(map)
    ax = gca;
    ax.XTickLabel = '';
    ax.YTickLabel = '';
    %ax.CLim = [0 1];
    set(gca, 'XColor', 'none')
    set(gca, 'YColor', 'none')
    grid on;
    ax.XTick = 0.5:5:39.5;
    ax.YTick = 0.5:5:39.5;
    set(gca, 'PlotBoxAspectRatio', [1 1 1])
    for j = ax.XTick
        line([j j], [0 40.5], 'Color', 'White', 'LineWidth', 1)
        line([0 40.5], [j j], 'Color', 'White', 'LineWidth', 1)
    end
  
end
clear ax j latency sub sub_rdm t
%% Calculate the number of trials per each condition in MEG dataset

n_trials = zeros(n,n_conditions);
DataDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/proc_data';

SJs = {'AA', 'AL', 'AR', 'BJB', 'CW', 'DJ', 'EC', 'FSM'...
    'JA', 'JC', 'JP', 'JS', 'LS', 'MC', 'NA', 'NC', 'SL', ...
    'SM', 'TK', 'TL','AC', 'AW', 'NM', 'SF'};

fname = 'realExemplarERF_bc.mat';

for sub = 1:n
    load(fullfile(DataDir, SJs{sub}, fname));
    for condition = 1:n_conditions
        if ~isnan(realExemplarERF{1,condition})
            n_trials(sub,condition) = size(realExemplarERF{1,condition},1);
        else
        end
    end
end

%save('/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/eLife_Wu2022/Figures/Data/n_trials_per_image.mat', 'n_trials')
clear condition DataDir fname realExemplarERF SJs sub 
%% Correlation between trial number and dissimilarity at the group-averaged level
% at a given timepoints

sub_rdm = nan(n_conditions,n_conditions,n);
latency = 0.47; %[0 0.2, 0.29, 0.47] ;
t = find(abs(times-latency(1)) < 0.0001);

for sub=1:n
    sub_rdm(:,:,sub) = squareform(Subjects(t,:,sub));
    for c = 1:n_conditions
        sub_rdm(c,c,sub) = nan;
    end
end

group_rdm = nanmean(sub_rdm,3); 
mean_distances = nanmean(group_rdm,1); % mean distances between condition n and all other conditions
mean_trialnumber = mean(n_trials); % mean trial number for each condition
SE_trialnumber = std(n_trials) ./ sqrt(n); % standard error of mean trial number

%mean_distances([25, 26]) = [];
%mean_trialnumber([25,26]) = [];

figure
scatter(mean_trialnumber, mean_distances, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b'); %red
hold on
h1=lsline;
h1.LineWidth=2;
h1.Color =  [1 0.4980 0.0549];
set(gca, 'PlotBoxAspectRatio', [1 1 1])
xlabel('Mean trial number per condition', 'FontSize', 10, 'Fontweight', 'normal')
ylabel('Mean dissimilarity to other conditions', 'FontSize', 10, 'Fontweight', 'normal')

%[r, pval]= corrcoef(mean_distances, mean_trialnumber)

figure(round(100*rand(1)))
hold on

for i = 1:length(mean_trialnumber)
    h=bar(i,mean_trialnumber(i));
    if ismember(i,1:20)
        set(h,'FaceColor',[0.4660 0.6740 0.1880]);
    elseif ismember(i,21:40)
        set(h,'FaceColor', [0.6350 0.0780 0.1840]);
    end
end
ax = gca;
ax.XLim = [0.5 40.5];
ax.YLim = [0 16];
%xlabel('condition number')
ylabel('Number of trials', 'FontSize', 10, 'Fontweight', 'normal')
pbaspect([3.5,1,1])
errorbar(1:n_conditions, mean_trialnumber, SE_trialnumber, 'LineStyle', 'none', 'Color', 'k')
box off
clear ax c group_rdm h latency mean_distances mean_trialnumber pval r SE_trialnumber sub sub_rdm t 
%%
latency = 0.47; %[0 0.2, 0.29, 0.47] ;
t = find(abs(times-latency(1)) < 0.0001);

% Subjects' RDM at a given time point 
sub_rdm = nan(40,40,24);
for s=1:24
    sub_rdm(:,:,s) = squareform(Subjects(t,:,s));
    for c = 1:40
        sub_rdm(c,c,s) = nan;
    end
end

%group_rdm = nanmean(sub_rdm,3); 
%mean_distances = nanmean(group_rdm,1);
%mean_distances([25, 26]) = [];
%mean_trialnumber = mean(n_trials);
%mean_trialnumber([25,26]) = [];
R=nan(24,1);
P=nan(24,1);
figure
for i = 1:24
    subplot(4,6,i)
    
    this_rdm = sub_rdm(:,:,i);
    mean_distances = nanmean(this_rdm,1);
    this_trialnumber = n_trials(i,:);
    
    missed = find(isnan(mean_distances));
    mean_distances(missed) = [];
    this_trialnumber(missed) = [];
     
    h = scatter(this_trialnumber, mean_distances, 'o', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'b'); %red
    hold on
    h1=lsline;
    h1.LineWidth=2;
    h1.Color =  [1 0.4980 0.0549];
    set(gca, 'PlotBoxAspectRatio', [1 1 1])
    
    [r, pval]= corrcoef(mean_distances, this_trialnumber);
    R(i,1)=r(2,1); P(i,1)=pval(2,1);
    clear h h1 this_rdm mean_distances this_trilannumber r pval missed
    %xlabel('Mean trial number per condition', 'FontSize', 10, 'Fontweight', 'normal')
    %ylabel('Mean dissimilarity to other conditions', 'FontSize', 10, 'Fontweight', 'normal')
end
[p, h, stats] = signrank(R, 0, 'Tail', 'both')
figure; histfit(R,10)

%%

R=nan(1,251);
Pval=nan(1,251);

for t= 1:251
    sub_rdm = nan(40,40,24);
    for s=1:24
        sub_rdm(:,:,s) = squareform(Subjects(t,:,s));
        for c = 1:40
            sub_rdm(c,c,s) = nan;
        end
    end
    
    group_rdm = nanmean(sub_rdm,3);
    mean_distances = nanmean(group_rdm,1);
    mean_distances([25, 26]) = [];
    mean_trialnumber = mean(n_trials);
    mean_trialnumber([25,26]) = [];
    
    
    [r, p]= corrcoef(mean_distances, mean_trialnumber);
    R(t)=r(2,1);Pval(t)=p(2,1);
    clear sub_rdm group_rdm mean_distances mean_trialnumber r p
end

figure;
plot(times, R)
pbaspect([2 1.3 1])
hold on
   
if any(Pval<0.05)
        sig_time = times(Pval<0.05);
        plot(sig_time, 0.4, 'Marker', 's',...
            'Markersize', 3, 'MarkerFaceColor', 'blue', 'MarkerEdge', 'none');
end
ax = gca;
ax.XLim = [-0.5 2];
ax.YLim = [-0.8 0.5];
ax.YTick = -0.8:0.2:0.4;
plot(times, zeros(1, length(times)), 'LineStyle', '--', 'color', 'k');
%line([sig_time(1) sig_time(1)],ax.YLim, 'LineStyle', '--', 'color', 'k')
%ax.YTickLabel = {'46','48', '50', '52', '54', '56'};
line([0 0], ax.YLim, 'color', 'k')

xlabel('Time (sec) relative to stimulus onset', 'FontSize', 10, 'Fontweight', 'normal')
ylabel('#trial - dissimilarity correlation ', 'FontSize', 10, 'Fontweight', 'normal')
box off
clear ax c Pval R s Sig_time t
%% Plot the across-subject 
R=nan(24, 251);
P=nan(24, 251);

%mean_distances=nan(24, 251);
%sub_trialnumber=nan(24, 251);

for t=1:251
    for s=1:24
        sub_rdm = squareform(Subjects(t,:,s));
        for c = 1:40
            sub_rdm(c,c) = nan;
        end
        mean_distances = nanmean(sub_rdm,1);
        sub_trialnumber = n_trials(s,:);
        missed = find(isnan(mean_distances));
        mean_distances(missed) = [];
        sub_trialnumber(missed) = [];
        [r, pval]= corrcoef(mean_distances, sub_trialnumber);
        R(s,t)=r(2,1); P(s,t)=pval(2,1);
    end  
end

H = nan(1,251);
Pval=nan(1,251);

for i_time = 1:length(times)   
    [p,h,~] = signrank(R(:,i_time), 0, 'tail','both');
    Pval(i_time) = p; H(i_time) = h;
    clear p h stats
end

figure;
mu = nanmean(R,1);
SE = nanstd(R)./sqrt(24);
shadedErrorBar(times, mu, SE, 'lineprops', {'color', 'b', 'LineWidth', 2});      
pbaspect([2 1.3 1])
hold on
   
if any(Pval<0.05)
        sig_time = times(Pval<0.05);
        plot(sig_time, 0.4, 'Marker', 's',...
            'Markersize', 3, 'MarkerFaceColor', 'blue', 'MarkerEdge', 'none');
end

ax = gca;
ax.XLim = [-0.5 2];
ax.YLim = [-0.5 0.5];
ax.YTick = -0.4:0.2:0.4;
plot(times, zeros(1, length(times)), 'LineStyle', '--', 'color', 'k');
xlabel('Time (sec) relative to stimulus onset', 'FontSize', 10, 'Fontweight', 'normal')
ylabel('#trial - dissimilarity correlation ', 'FontSize', 10, 'Fontweight', 'normal')
line([0 0], ax.YLim, 'color', 'k')
clear pval Pval r R s SE sig_time sub_rdm sub_trialnumber ax c H i i_time mean_distances y1 t missed mu P