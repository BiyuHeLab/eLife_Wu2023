% event_dict = {'face1_r', 'face2_r', 'face3_r', 'face4_r', 'face5_r', 'face6_r',...
%               'house1_r', 'house2_r', 'house3_r', 'house4_r', 'house5_r', 'house6_r',...
%               'object1_r', 'object2_r', 'object3_r', 'object4_r', 'object5_r', 'object6_r',...
%               'animal1_r', 'animal2_r', 'animal3_r', 'animal4_r', 'animal5_r', 'animal6_r',...
%               'face1_nr', 'face2_nr', 'face3_nr', 'face4_nr', 'face5_nr', 'face6_nr',...
%               'house1_nr', 'house2_nr', 'house3_nr', 'house4_nr', 'house5_nr', 'house6_nr',...
%               'object1_nr', 'object2_nr', 'object3_nr', 'object4_nr', 'object5_nr', 'object6_nr',...
%               'animal1_nr', 'animal2_nr', 'animal3_nr', 'animal4_nr', 'animal5_nr', 'animal6_nr',...
%               'noresponse'};
function [realExemplarERF, scrExemplarERF] = organize_ERF(datain)
%% organize MEG time-series data based on exemplars and recognition reports

nConditions = 49;
scrambled = 6:6:48;
motor = 49;
real= 1:nConditions;
real([scrambled motor])= [];
ERF = cell(40, nConditions);

for t = 1:length(datain.time)
    for exemplar = 1:nConditions
        tr_num = find(datain.trialinfo==exemplar);
        if ~isempty(tr_num)
            tmp = datain.trial(tr_num, :,t);
        else
            tmp = NaN;
        end
        ERF{t, exemplar}=tmp;
        clear tmp tr_num
    end
end
realExemplarERF = ERF(:, real);
scrExemplarERF = ERF(:, scrambled);

