
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

elseif strcmpi(ModelName, 'Animacy_1')
    ModelRDM = [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0];
    ModelRDM = kron(ModelRDM, ones(5,5));
    ModelRDM = [ModelRDM, ones(20,20); ones(20,20), ModelRDM];
    for i = 1:40
        ModelRDM(i,i) = 0;
    end
    
elseif strcmpi(ModelName, 'Category_1') 
    ModelRDM = ones(20);
    for i = [1 6 11 16]
        ModelRDM(i:i+4, i:i+4) = 0;
    end  
    ModelRDM = [ModelRDM, ones(20); ones(20), ModelRDM];
    for i= 1:40
        ModelRDM(i,i) = 0;
    end      
 
elseif strcmpi(ModelName, 'Rec_Category_1')
    RecQ = ones(20,20);
    
    for i = [1 6 11 16]
        RecQ(i:i+4, i:i+4) = 0;
    end
    ModelRDM = [RecQ, ones(20,20); ones(20,20), zeros(20,20)];
    for i = 1:40
        ModelRDM(i,i) = 0;
    end
    clear RecQ

elseif strcmpi(ModelName, 'Rec_Animacy_1')
    ModelRDM = [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0];
    ModelRDM = kron(ModelRDM, ones(5,5));
    ModelRDM = [ModelRDM, ones(20,20); ones(20,20), zeros(20,20)];
    for i = 1:40
        ModelRDM(i,i) = 0;
    end
    
elseif strcmpi(ModelName, 'Unrec_Animacy_1')
    ModelRDM = [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0];
    ModelRDM = kron(ModelRDM, ones(5,5));
    ModelRDM = [zeros(20,20), ones(20,20); ones(20,20), ModelRDM];
    for i = 1:40
        ModelRDM(i,i) = 0;
    end
       
elseif strcmpi(ModelName, 'Unrec_Category_1')
    RecQ = ones(20,20);
    
    for i = [1 6 11 16]
        RecQ(i:i+4, i:i+4) = 0;
    end
    ModelRDM = [zeros(20,20), ones(20,20); ones(20,20), RecQ];
    for i = 1:40
        ModelRDM(i,i) = 0;
    end
    clear RecQ
    
elseif strcmpi(ModelName, 'Recognition')
       ModelRDM = zeros(20,20);
       ModelRDM = [ModelRDM, ones(20,20); ones(20,20), ones(20,20)];
       for i= 1:40
            ModelRDM(i,i) = 0;
       end    

elseif strcmpi(ModelName, 'Rec_Animacy_2')
    ModelRDM = [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0];
    ModelRDM = kron(ModelRDM, ones(5,5));
    ModelRDM = [ModelRDM, ones(20,20); ones(20,20), ones(20,20)];
    for i = 1:40
        ModelRDM(i,i) = 0;
    end    
       
 elseif strcmpi(ModelName, 'Rec_Category_2')
    RecQ = ones(20,20);
    
    for i = [1 6 11 16]
        RecQ(i:i+4, i:i+4) = 0;
    end
    ModelRDM = [RecQ, ones(20,20); ones(20,20), ones(20,20)];
    for i = 1:40
        ModelRDM(i,i) = 0;
    end
    clear RecQ      

elseif strcmpi(ModelName, 'Unrec_Clustering')
       ModelRDM = zeros(20,20);
       ModelRDM = [ones(20,20), ones(20,20); ones(20,20), ModelRDM];
       for i= 1:40
            ModelRDM(i,i) = 0;
       end

elseif strcmpi(ModelName, 'Unrec_Animacy_2')
    ModelRDM = [0 1 1 0; 1 0 0 1; 1 0 0 1; 0 1 1 0];
    ModelRDM = kron(ModelRDM, ones(5,5));
    ModelRDM = [ones(20,20), ones(20,20); ones(20,20), ModelRDM];
    for i = 1:40
        ModelRDM(i,i) = 0;
    end  
       
  elseif strcmpi(ModelName, 'Unrec_Category_2')
    RecQ = ones(20,20);
    
    for i = [1 6 11 16]
        RecQ(i:i+4, i:i+4) = 0;
    end
    ModelRDM = [ones(20,20), ones(20,20); ones(20,20), RecQ];
    for i = 1:40
        ModelRDM(i,i) = 0;
    end
    clear ReQ 
      
elseif strcmpi(ModelName, 'Category_2') 
    ModelRDM = ones(20);
    for i = [1 6 11 16]
        ModelRDM(i:i+4, i:i+4) = 0;
    end  
    ModelRDM = repmat(ModelRDM,2,2);
    for i = 1:40
        ModelRDM(i,i) = 0;
    end
end

ModelVector = squareform(ModelRDM);
ModelVector = ModelVector';


