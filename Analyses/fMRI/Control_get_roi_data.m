function [beta_hat, RoiSize] = Control_get_roi_data(SJ, SubEVs, UserOptions, MaskName)

ProjDir     = UserOptions.ProjDir;
%RoiType     = UserOptions.RoiType;
if strcmp(UserOptions.RoiSpace, 'native')
    RoiSpace    = 'native';
elseif strcmp(UserOptions.RoiSpace, '_standard')
    RoiSpace    = 'standard';
end
%%
RoiMasks = {};
RoiVoxels = [];

subj_dir = fullfile(ProjDir, 'ROI_masks', ['sub' SJ]);
RoiMasks = fullfile(subj_dir,[MaskName '_native.nii.gz']);


% if ROI is not existent, all remaining statements in the body of the
% loop will be skipped and continue with the next iteration
if isempty(RoiMasks)
    RoiSize = nan;
    beta_hat = nan;
    u_hat = nan;
    Sw_hat = nan;
    return
else
       
    if isstruct(RoiMasks)
        flds = fields(RoiMasks);
        for i = 1:numel(flds)
            roi = load_nifti(RoiMasks.(flds{i}));
            RoiVoxels = [RoiVoxels find(roi.vol)'];
        end
    elseif ischar(RoiMasks)
        roi = load_nifti(RoiMasks);
        RoiVoxels = find(roi.vol)';
    end
    %make sure that there are no overlapping voxels.
    %Values in the output matrix represent the voxel indices
    RoiVoxels = unique(RoiVoxels); 
end

clear RoiMasks RoiName RoiType roi roi_dir RoiSpace MaskName fleds i

%%
for ev = 1:size(SubEVs,1)
    % identify blocks wherein the given EV was modelled, starting with
    % EV1 which is face1_recognized
    blocks = find(SubEVs(ev,:));
    % create a matrix with the dimension nSamples x nVoxels
    beta_hat{1, ev} = zeros(length(blocks), length(RoiVoxels));
    
    
    for i = 1:length(blocks)
        % find the EPI and extract the beta RoiVoxels = unique(RoiVoxels);estimates from voxels within
        % the ROI
        EPIDir = '/isilon/LFMI/VMdrive/data/HLTP_fMRI';
        EPI = fullfile(EPIDir, ['sub' SJ], 'proc_data/func', ['block' num2str(blocks(i))], 'StimIdentityGLM.feat/stats', ['cope' num2str(ev) '_2highres.nii']);
        tmp = load_nifti(EPI);
        beta_hat{1, ev}(i,:) = tmp.vol(RoiVoxels);
        clear EPI tmp
    end
    clear blocks
end
RoiSize = length(RoiVoxels);

