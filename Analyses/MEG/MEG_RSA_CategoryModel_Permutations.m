function [PermRho, PermStats] = MEG_RSA_CategoryModel_Permutations(ModelVector, SubRDMs, n_perms)
% This function permutes exemplar labels in recognized and unrecognized
% trials, respectively.
% Then, it computes Spearman's rhos between the permuted data and a model
% RDM and carries out statistical testing for each permutation
% Last update: 3/7/2023 @Yuan-hao Wu

for i_perm = 1:n_perms
    tic
    rp = (randperm(size(SubRDMs.seen,1)))'; % use the same randomization across all ROIs in each permuation step
    
    SeenPerm = SubRDMs.seen(rp,:,:);
    UnseenPerm = SubRDMs.unseen(rp,:,:);
    
    % compute time-varying data-model correlations
    % separately
    for sub = 1:size(SubRDMs.seen,3)
        [r, ~] = corr(SeenPerm(:,:,sub), ModelVector, 'type', 'Spearman', 'rows', 'complete', 'tail', 'right');
        PermRho.seen(sub,:, i_perm)= r;
        clear r
        
        [r, ~] = corr(UnseenPerm(:,:,sub), ModelVector, 'type', 'Spearman', 'rows', 'complete', 'tail', 'right');
        PermRho.unseen(sub,:, i_perm)= r;
        clear r
    end
    clear SeenPerm UnseenPerm
    
    
    for t = 1:size(SubRDMs.seen,2)
        [p,h,stats] = signrank(PermRho.seen(:,t,i_perm), 0, 'tail','right');
        PermStats.seen.p(i_perm,t) = p; PermStats.seen.h(i_perm,t) = h;
        PermStats.seen.zval(i_perm,t) = stats.zval; PermStats.seen.signedrank(i_perm,t) = stats.signedrank;
        clear p h stats
        
        [p,h,stats] = signrank(PermRho.unseen(:,t, i_perm), 0, 'tail','right');
        PermStats.unseen.p(i_perm,t) = p; PermStats.unseen.h(i_perm,t) = h;
        PermStats.unseen.zval(i_perm,t) = stats.zval; PermStats.unseen.signedrank(i_perm,t) = stats.signedrank;
        clear p h stats
    end
    
    elapsed_time  = toc;
    display(['Permutation iteration ' num2str(i_perm, '%04.f') ' completed, Elapsed time: ' num2str(elapsed_time)])
    clear rp t elapsed_time
end



