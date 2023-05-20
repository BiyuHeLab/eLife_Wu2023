clear
SJs = {'AA', 'AC','AL', 'AR', 'AW', 'BJB', 'CW', 'DJ', 'EC', 'FSM', ...
    'JA', 'JC', 'JP', 'JS', 'LS', 'MC', 'NA', 'NM', 'NC', 'SF','SL', ...
    'SM', 'TK', 'TL'};

DataDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/distance_data';
DistanceMeasure = 'realAvg1-Pearson_bc_400Hz';
times = -0.5:0.0025:2;
nConditions = 40;
nDistances = (nConditions * nConditions - nConditions) / 2;


%%
tmp = nan(length(times), nDistances, length(SJs));

for s = 1:length(SJs)
    load(fullfile(DataDir, SJs{s}, 'Unnormalized_Data', [DistanceMeasure '.mat']));
    tmp(:,:,s) = AvgDistances;
    clear AvgDistances
end

AvgData = nanmean(tmp,3);
save(fullfile(DataDir, 'Grand_Average', 'Unnormalized_Data', [DistanceMeasure '.mat']), 'AvgData')