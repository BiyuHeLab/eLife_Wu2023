function Dat = import_mne_epochs(SubDir, EpoName)
%This function imports epochs in mne to fieldtrip format. It requires
%-epo.fif and -eve.fif being stored in the same folder. The first input
%argument is the path for the stored mne files, whereas the second input is
%the file name without the mne extension.


fiff_file = fullfile(SubDir, [EpoName '-epo.fif']);
events_file= fullfile(SubDir, [EpoName '-eve.fif']);

cfg = [];
cfg.dataset = fiff_file;
Dat = ft_preprocessing(cfg);

event_file = mne_read_events(events_file);
Dat.trialinfo = event_file(:,3);
Dat.cfg.trl(:,4) = event_file(:,3);

