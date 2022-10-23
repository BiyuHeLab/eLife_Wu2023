%==========================================================================
%This script imports mne epochs to fieldtrip format, cuts the imported 
%epochs to desired length and applies baseline correction to them.
%In addition, the epochs will be re-organized into 40 conditions for the
%later RSA.
% requires fieldtrip toolbox, 'import_mne_epochs.m', and 'organize_ERF.m'
%==========================================================================
clear;clc
addpath('/isilon/LFMI/VMdrive/YuanHao/toolboxes/fieldtrip-master')
addpath('/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/Scripts/DataPreparation')
ft_defaults 

SJs = {'AA', 'AC', 'AL', 'AR', 'AW', 'BJB', 'CW', 'DJ', 'EC', 'FSM'...
   'JA', 'JC', 'JP', 'JS', 'LS', 'MC', 'NA', 'NC', 'NM', 'SF', 'SL', ...
   'SM', 'TK', 'TL'};

DataDir='/isilon/LFMI/VMdrive/YuanHao/HLTP_Fusion/MEG/proc_data';
EpoName= 'HLTP_0-35Hz_400Hz';
%%
for sj = 1:length(SJs)
    SubDir = fullfile(DataDir, SJs{sj});
    %import epochs from mne and save it in fieldtrip format
    Dat = import_mne_epochs(SubDir ,EpoName);
    save(fullfile(SubDir, [EpoName '.mat']), 'Dat')
    clear Dat
% ===================================================================
% shorten epochs and apply baseline correction
% ===================================================================
    load(fullfile(SubDir, [EpoName '.mat']));
    cfg = [];
    cfg.toilim = [-0.5 2]; %set the start and the end of a trial
    tmp = ft_redefinetrial(cfg, Dat);
    
    cfg= [];
    cfg.baseline = [-0.5 0]; %'no'% set the start and the end of the baseline
    Dat_bc = ft_timelockbaseline(cfg, tmp);
    save(fullfile(SubDir, [EpoName '.mat']), 'Dat')
    clear tmp Dat cfg Dat_bc

% ================================================== 
% Reorganize the data for RSA 
% ==================================================    
    load(fullfile(SubDir, [EpoName '.mat']));
    [realExemplarERF, scrExemplarERF] = organize_EPF(Dat_nobc);
    allExemplarERF = [realExemplarERF, scrExemplarERF];
    
    save(fullfile(SubDir, 'allExemplarERF_nobc.mat'), 'allExemplarERF');
    save(fullfile(SubDir, 'realExemplarERF_bc.mat'), 'realExemplarERF');
    save(fullfile(SubDir, 'scrExemplarERF_bc.mat'), 'scrExemplarERF');
        
    clear Dat_bc allExemplarERF realExemplarERF scrExemplarERF SubDir 
end