% Prepare data for decoding analysis
% transform the ExemplarERF to a format optimized for decoding analysis
clear;clc
MEG_DataDir = '/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/proc_data/';

SJs = {'AA', 'AL', 'AR', 'BJB', 'CW', 'DJ', 'EC', 'FSM'...
    'JA', 'JC', 'JP', 'JS', 'LS', 'MC', 'NA', 'NC', 'SL', ...
    'SM', 'TK', 'TL','AC', 'AW', 'NM', 'SF'};

i_Recog = [1:5; 6:10; 11:15; 16:20];
i_Unrecog = i_Recog + 20;
Categories = {'Face', 'House', 'Object', 'Animal'};

for s = 1:length(SJs)
    tic
    FileName =  fullfile(MEG_DataDir, SJs{s}, 'realExemplarERF_bc_400Hz.mat');
    Data = load(FileName); Data = Data.realExemplarERF;
    
    %find conditions with 0 trial. Replace NAN with [];
    nantrials = find(cellfun(@(C) isequaln(C, NaN), Data));
    for i = 1:length(nantrials)
        Data{nantrials(i)} = [];
    end
     
    for i_time = 1:size(Data,1)
        for i = 1:length(Categories)        
            RecogCat.(Categories{i})(:,:,i_time) = ...
                [Data{i_time,i_Recog(i,1)};
                Data{i_time,i_Recog(i,2)};
                Data{i_time,i_Recog(i,3)};
                Data{i_time,i_Recog(i,4)};
                Data{i_time,i_Recog(i,5)}];
            
            UnrecogCat.(Categories{i})(:,:,i_time) = ...
                [Data{i_time,i_Unrecog(i,1)};
                Data{i_time,i_Unrecog(i,2)};
                Data{i_time,i_Unrecog(i,3)};
                Data{i_time,i_Unrecog(i,4)};
                Data{i_time,i_Unrecog(i,5)}];           
        end
    end
    save(fullfile(MEG_DataDir, SJs{s}, 'DataForDecoding.mat'), 'RecogCat', 'UnrecogCat');
    clear FileName Data RecogCat UnrecogCat nantrials
    elapsed_time = toc;
    disp(['Sj ' num2str(s, '%02.f') ' : ' num2str(elapsed_time) ' _sec']);
end


