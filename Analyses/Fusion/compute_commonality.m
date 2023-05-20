
function commModel  = compute_commonality(megRDM,fmriRDM, modelRDM)
 % Estimate commonality coefficient between three variables (MEG RDM from a
 % given time point, fMRI RDM from a given location, and a given conceptual
 % model).
 % The standard equation (Seibold & McPhee, 1979) is: 
 % Rsquared(MEG.fMRI) + Rsquared(MEG.Model) - Rsquared(MEG.fMRI, Model)
 % Here we use a computationally more efficient approach (as in
 % Hebart et al., 2018 eLife), which is
 % bivariate correlation(MEG, MRI) - semipartialcorrelation(MEG,MRI.model)
 % Author: @Yuan-hao Wu
 % Last update: 9/25/2022

commModel = zeros(1,size(megRDM,2));

for i_time = 1:size(megRDM,2)
    thisMegRDM = megRDM(:, i_time);
    rMRI_MEG = correlate([thisMegRDM fmriRDM],'type','spearman','method','corr');
    rMRI_MEG_model = correlate([thisMegRDM fmriRDM modelRDM],'type','spearman','method','semipartialcorr');
    commModel(i_time) = rMRI_MEG(1,2)^2 - rMRI_MEG_model(1,2)^2;
    
    %below the standard implementation a la Seibold and McPhee  
    %Rsquared12 = (corr(thisMegRDM, fmriRDM,'type','spearman'))^2;
    %Rsquared13 = (corr(thisMegRDM, modelRDM,'type','spearman'))^2;
    %Rsquared123 = multiRegress2var(thisMegRDM, fmriRDM, modelRDM);
    %commModel(i_time) = Rsquared12 + Rsquared13  - Rsquared123;   
end

