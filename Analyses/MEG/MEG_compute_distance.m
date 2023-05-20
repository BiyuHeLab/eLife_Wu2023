% Calculate the 1-Pearson's r distances between the averaged ERF patterns
% elicited by each exemplar separately
% 1.) calculate the mean ERF pattern for each condition.
% 2.) compute the 1-Pearson's r distance between all 
% condition pairs: 40 conditions yield 780 pair-wise distances.
% n_distances = (Nconditions^2 - Nconditions)/2
% output n_timepoints x n_distances matrix for each subject

clear; clc;
TrialType = 'real'; %['real', 'scr', 'all'];
BaselineCorrection = 'bc';
DistanceMeasure    = 'correlation'; %e.g., 'correlation', 'euclidean', 'seuclidean',...
EPFDir            = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/proc_data'; % location of ERF patterns
OutputDir      = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/distance_data';
NoiseNormalization  = 'none';
fname = [TrialType 'ExemplarERF_' BaselineCorrection '_400Hz'];
times = -0.5:0.0025:2;
nTimepoints = length(times);

if strcmp(TrialType, 'all') ==1
   conds = 1:48;
elseif strcmp(TrialType, 'real') ==1
    conds = 1:40;
elseif strcmp(TrialType, 'scr') ==1
    conds = 41:48;
end
nConditions = length(conds);
%nChannels   = 272
nDistances  = (nConditions*nConditions-nConditions)/2;
%************************************************************************
SJs = {'AA', 'AC', 'AL', 'AR', 'AW', 'BJB', 'CW', 'DJ', 'EC', 'FSM',...
       'JA', 'JC', 'JP', 'JS', 'LS', 'MC', 'NA', 'NC', 'NM', 'SF', 'SL', ...
       'SM', 'TK', 'TL'};
%%
for subj = 1:length(SJs)
    %load subject-specific ERF patterns
    filepath = fullfile(EPFDir, SJs{subj}, [fname '.mat']);
    ExemplarERF = load(filepath); ExemplarERF = ExemplarERF.realExemplarERF;
    
    nChannels = size(ExemplarERF{1,nConditions},2);   
    AvgDistances = zeros(nTimepoints, nDistances);   
    
    % Loop over all timepoints
    for i_time = 1:nTimepoints
        AvgERF = zeros(nConditions, nChannels);
        t_data = ExemplarERF(i_time,conds);     % determine the data of this particualr timepoint
        
        % compute the mean ERF activity pattern of each exemplar, respectively 
        for i = 1:nConditions
            AvgERF(i,:) = nanmean(t_data{i},1);  
        end
        % now compute the between-exemplar dissimilarities
        AvgDistances(i_time, :) = pdist(AvgERF, DistanceMeasure);
        clear t_data AvgERF
    end
    save(fullfile(OutputDir, SJs{subj}, 'Unnormalized_Data', [TrialType 'Avg1-Pearson_' BaselineCorrection '.mat']), 'AvgDistances')
    clear SubjectData ExemplarERF AvgDistances nChannels i_time filepath
end                        