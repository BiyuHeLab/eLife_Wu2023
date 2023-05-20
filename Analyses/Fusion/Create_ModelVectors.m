
function ModelVector = Create_ModelVectors(ModelName)
% Create a n x 1 vector corresponding to a conceptual model RDM
% n is the number of condition pairs and derived from
% (nConditions x nConditions - nConditions) / 2
% Author: @Yuan-hao Wu
% Last Update: 9/25/2022

if strcmpi(ModelName, 'TwoState') 
    RecStates = zeros(20,20);
    Mixed = ones(20,20);
    ModelRDM = [RecStates, Mixed; Mixed, RecStates];
    for i= 1:40
            ModelRDM(i,i) = 0;
       end
    clear RecStates Mixed

  
elseif strcmpi(ModelName, 'Recognition')
       ModelRDM = zeros(20,20);
       ModelRDM = [ModelRDM, ones(20,20); ones(20,20), ones(20,20)];
       for i= 1:40
            ModelRDM(i,i) = 0;
       end    
end

ModelVector = squareform(ModelRDM);
ModelVector = ModelVector';
