clear; clc
method = 'ensemble_trial_5folds';

addpath('/isilon/LFMI/VMdrive/YuanHao/toolboxes/decoding_toolbox_v3.997')

MEG_DataDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/proc_data/';
DecodingDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/decoding_data/';

SJs = {'AA', 'AC', 'AL', 'AR', 'AW,', 'BJB', 'CW', 'DJ', 'EC', 'FSM'...
    'JA', 'JC', 'JP', 'JS', 'LS', 'MC', 'NA', 'NC', 'NM', 'SL', ...
    'SF', 'SM', 'TK', 'TL'};

class = {'seen', 'unseen'}; %'House'; , 'Animal'};
obj = {'Face', 'House', 'Object', 'Animal'};
times = -0.5:0.01:2;
%%
for s =1:length(SJs)
    FileName =  fullfile(MEG_DataDir, SJs{s}, 'DataForDecoding.mat');
    Data = load(FileName);
    SeenObj = Data.RecogCat; UnseenObj = Data.UnrecogCat;
    clear Data
    
    data1 = [];
    data2 = [];
    for i = 1:length(obj)
        data1 = [data1; SeenObj.(obj{i})];
        data2 = [data2; UnseenObj.(obj{i})];
    end
    clear SeenObj UnseenObj
    %% Randomly assign trials to 5 cv folds without replacement.
    nfolds = 5;
    labels = [];
    chunks = [];
    idx_perm = [];
    ntrials = [size(data1,1), size(data2,1)];
    
    for i = 1:size(class,2)
        % determine how many trials goe to each single cv fold
        nRemainders = rem(ntrials(i),nfolds);
        nt=(ntrials(i)-nRemainders)/nfolds;
        nt= repmat(nt, nfolds, 1);
        if nRemainders ~=0
            nt(1:nRemainders,1) = nt(1:nRemainders,1) +1;
        end       
        % create vectors to indicate each trial's label and chunk
        for j=1:nfolds
            chunks = [chunks; j*ones(1, nt(j))'];
        end
        labels= [labels; repmat(i , ntrials(i), 1)];
        
        %generate randomized trial order for each object category
        idx_perm{i} =randperm(ntrials(i))';
        clear nRemainders nt j
    end
    
    %radomized trial order
    data1 = data1(idx_perm{1},:,:);
    data2 = data2(idx_perm{2},:,:);
    clear i idx_perm ntrials nfolds
    %% parameter settings for TDT
    % Set defaults for TDT
    cfg = decoding_defaults;
    cfg.files.mask = '';
    cfg.analysis = 'wholebrain';
    cfg.files.label = labels;
    cfg.files.chunk = chunks;
    
    % save a description
    all_chunks = unique(cfg.files.chunk);
    all_labels = unique(cfg.files.label);
    
    ct = zeros(length(all_labels),length(all_chunks));
    for ifile = 1:length(cfg.files.label)
        curr_label = cfg.files.label(ifile);
        curr_chunk = cfg.files.chunk(ifile);
        f1 = all_labels==curr_label; f2 = all_chunks==curr_chunk;
        ct(f1,f2) = ct(f1,f2)+1;
        cfg.files.name(ifile,1) = {sprintf('class%i_block%i_%i', curr_label, curr_chunk, ct(f1,f2))};
    end
    cfg.scale.method = 'z';
    %cfg.scale.estimation = 'all';
    cfg.scale.check_datatrans_ok = true;
    %cfg.feature_transformation.method = 'PCA';
    %cfg.feature_transformation.estimation = 'all';
    %cfg.feature_transformation.critical_value = 0.01;
     
    cfg.design = make_design_cv(cfg);
    cfg.plot_design = 1;
    display_design(cfg);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Balance ensemble approach: Subsample the more frequent class
    % repeatedly, run multiple classification iterations and predict all labels
    % for each of classifiers. Use combined decision values to create a
    % majority vote of all classifiers for one final prediction
    % PRO: uses all data, better performance than repeated subsampling
    % CON: may become slow
    % Implementation:
    cfg.decoding.method = 'classification';
    cfg.decoding.train.classification.model_parameters = '-s 0 -t 2 -c 1 -b 0 -q';
    cfg.design.unbalanced_data = 'ok';
    cfg.decoding.software = 'ensemble_balance';
    cfgd = decoding_defaults; % to use default values
    cfg.decoding.train.classification.model_parameters = [];
    cfg.decoding.train.classification.model_parameters.software = 'libsvm';
    cfg.decoding.train.classification.model_parameters.n_iter = 100;
    cfg.decoding.train.classification.model_parameters.model_parameters = cfgd.decoding.train.classification.model_parameters;
    cfg.decoding.test.classification.model_parameters = [];
    cfg.decoding.test.classification.model_parameters.model_parameters = cfgd.decoding.test.classification.model_parameters;
    cfg.results.output = {'accuracy_minus_chance', 'balanced_accuracy_minus_chance', ...
        'decision_values', 'predicted_labels', 'AUC_minus_chance', 'AUC_pairwise', 'dprime'};
    cfg.results.overwrite = 1;
    clear all_chunks all_labels chunks ct curr_chunk curr_label f1 f2 FileName i ifile j labels
    %% Run decoding for time-varying MEG data
    for i_time =1:length(times)
        
        data = [data1(:,:,i_time);
                data2(:,:,i_time)];
             
        % Set the output directory where data will be saved, e.g. 'c:\exp\results\buttonpress'
        cfg.results.dir = fullfile(DecodingDir, 'RecognitionState', SJs{s}, method, ['t' num2str(i_time, '%03.f')]);
        if exist(cfg.results.dir, 'dir') ~= 7
            mkdir(cfg.results.dir)
        else
        end
        cd(cfg.results.dir)
        
        % Prepare data for passing
        passed_data.data = data;
        passed_data.mask_index = 1:size(data, 2); % use all voxels
        passed_data.files = cfg.files;
        passed_data.hdr = ''; % we don't need a header, because we don't write img-files as output (but mat-files)
        passed_data.dim = [size(data, 2), 1, 1]; % add dimension information of the original data
        
        % Run decoding
        [results, cfg] = decoding(cfg, passed_data);
        clear data passed_data results
    end
    clear cfg cfgd data1 data2 idx_perm
end