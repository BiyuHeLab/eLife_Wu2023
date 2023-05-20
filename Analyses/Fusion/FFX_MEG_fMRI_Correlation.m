% computes the correlations between the 40x40 group-average fMRI RDM and
% the 40x40 group-average MEG RDM.

clear; clc
ProjDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/';
MRI_Dir = fullfile(ProjDir,'fMRI', 'ROI_Data', 'Grand_Average', 'Unnormalized_Data');
MEG_Dir = fullfile(ProjDir, 'MEG', 'distance_data', 'Grand_Average', 'Unnormalized_Data');

load(fullfile(MEG_Dir, 'realAvg1-Pearson_bc_400Hz.mat'), 'AvgData');
megRDM = AvgData; clear AvgData
megRDM = megRDM';
load(fullfile(MRI_Dir, 'realAvg1-Pearson_n23.mat'), 'AvgData');
ROIs = fieldnames(AvgData);
mriRDM = AvgData; clear AvgData

ROIs = ROIs(1:29);
%ROIs([9,10,22,23]) = [];


for r = 1:length(ROIs)
    tmp1 = zeros(1,1001);
    tmp2 = zeros(1,1001);
    roiRDM = mriRDM.(ROIs{r})';
    for i = 1:1001
        [tmp1(1,i), tmp2(1,i)] = corr(megRDM(:,i), roiRDM, 'type', 'Spearman',...
        'row', 'complete', 'tail', 'right');
    end
    AvgRho.(ROIs{r}) = tmp1;
    pval.(ROIs{r}) = tmp2;
    clear roi_data tmp1 tmp2
end

%save(fullfile(TargetDir, 'MEG-fMRI_Correlation_n20n23_bc_400Hz.mat'), 'AvgRho', 'pval');

times = -0.5:0.0025:2;

%% Randomization test with 5000 permutation samples
nPerms = 5000;
rho_perm = nan(nPerms, length(times),length(ROIs));

parfor r = 1:length(ROIs)
    tic
    for k = 1:nPerms
        roiRDM = fMRI_data.AvgData.(ROIs{r})';
        roiRDM = roiRDM(randperm(length(roiRDM)));
        rho_perm(k,:,r) = corr(megRDM,roiRDM, 'type', 'Spearman', 'row', 'complete');
        %clear fmriRDM
    end
    
    elapsed = toc;
    disp('                                                        ')
    disp('==========================================================================')
    disp([ROIs{r} ' done. Elpased time: ' num2str(elapsed) ' sec.'])
    disp('==========================================================================')
    disp('                                                        ')
    
end
clear modelRDM megRDM k m r elapsed
save('/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/Fusion/Data/commonality_400Hz/MEG-fMRI_Fusion_5000perms.mat', 'rho_perm')



times = -0.5:0.0025:2;
%%
colors = {[0, 0.4470, 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
figure(round(100*rand(1)))
count =1;
for r = 9:16 %length(ROIs)
    subplot(2,4,count)
    ax = gca;
    %AvgRho = nanmean(Spearman.(Rois{r}).rho,1);
    %SE = nanstd(Spearman.(Rois{r}).rho)./sqrt(length(SJs)-1);
    %shadedErrorBar(times,AvgRho, SE, 'lineprops', {'color', colors{3}, 'LineWidth', 2})
    plot(times, AvgRho.(ROIs{r}), 'LineWidth', 1)
    hold on
    pbaspect([1.5 1 1])
    ax=gca;
    ax.YLim = [-0.3 0.5];
    ax.XLim = [-0.5 2];
    plot(times, zeros(1, length(times)), 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--')
    line([0 0], ax.YLim, 'Color', 'k')  
    title(ROIs{r}, 'Interpreter', 'none')
    box off
    count = count+1;
    
end

%% First average RDMs across regions, then compute correlations with models
AvgData = [];
for r= 1:4
    roi_data = fMRI_data.AvgData.(ROIs{r})';
    AvgData = [AvgData roi_data];
end
AvgData = nanmean(AvgData,2);
EVC = corr(time_data, roi_data, 'type', 'Spearman', 'row', 'complete');
    
AvgData = [];
for r= 5:8
    roi_data = fMRI_data.AvgData.(ROIs{r})';
    AvgData = [AvgData roi_data];
end
AvgData = nanmean(AvgData,2);
Obj = corr(time_data, roi_data, 'type', 'Spearman', 'row', 'complete');

AvgData = [];
for r= 11:19
    roi_data = fMRI_data.AvgData.(ROIs{r})';
    AvgData = [AvgData roi_data];
end
AvgData = nanmean(AvgData,2);
TaskPos = corr(time_data, roi_data, 'type', 'Spearman', 'row', 'complete');

AvgData = [];
for r= 20:29
    roi_data = fMRI_data.AvgData.(ROIs{r})';
    AvgData = [AvgData roi_data];
end
AvgData = nanmean(AvgData,2);
DMN = corr(time_data, roi_data, 'type', 'Spearman', 'row', 'complete');
clear AvgData

RegionAvg = [EVC, Obj, TaskPos, DMN];


colors = {[0, 0.4470, 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0 0 0]};
figure(round(100*rand(1)))

for r = 1:4%length(Rois)
    %subplot(2,4,count)
    ax = gca;
    %AvgRho = nanmean(Spearman.(Rois{r}).rho,1);
    %SE = nanstd(Spearman.(Rois{r}).rho)./sqrt(length(SJs)-1);
    %shadedErrorBar(times,AvgRho, SE, 'lineprops', {'color', colors{3}, 'LineWidth', 2})
    plot(times, RegionAvg(:,r), 'Color', colors{r}, 'LineWidth', 2)
    hold on
end
    pbaspect([1.5 1 1])
    ax.YLim = [-0.15 0.4];
    plot(times, zeros(1, length(times)), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    line([0 0], ylim, 'Color', 'k')  
    title('fMRI-MEG Fusion')
    box off
    legend('EVC', 'Object-selective areas', 'Task positive network', 'DMN')
    ylabel('Spearman''s rho', 'Fontsize', 12)
xlabel('time', 'Fontsize', 12 )

%% Compute MEG-fMRI correlationfor each ROI separately, and then average across ROIs
Correlations = nan(251,29); 
for r = 1:29
    Correlations(:,r) = AvgRho.(ROIs{r}); 
end

EVC = nanmean(Correlations(:,1:4),2);
Obj = nanmean(Correlations(:,5:8),2);
TaskPos = nanmean(Correlations(:,11:19),2);
DMN = nanmean(Correlations(:,20:29),2);

RegionAvg = [EVC, Obj, TaskPos, DMN];

colors = {[0, 0.4470, 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0 0 0]};
figure(round(100*rand(1)))

for r = 1:4%length(Rois)
    %subplot(2,4,count)
    ax = gca;
    %AvgRho = nanmean(Spearman.(Rois{r}).rho,1);
    %SE = nanstd(Spearman.(Rois{r}).rho)./sqrt(length(SJs)-1);
    %shadedErrorBar(times,AvgRho, SE, 'lineprops', {'color', colors{3}, 'LineWidth', 2})
    plot(times, RegionAvg(:,r), 'Color', colors{r}, 'LineWidth', 2)
    hold on
end
    pbaspect([1.5 1 1])
    ax.YLim = [-0.15 0.4];
    plot(times, zeros(1, length(times)), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '--')
    line([0 0], ylim, 'Color', 'k')  
    title('fMRI-MEG Fusion')
    box off
    legend('EVC', 'Object-selective areas', 'Task positive network', 'DMN')
    ylabel('Spearman''s rho', 'Fontsize', 12)
xlabel('time', 'Fontsize', 12 )