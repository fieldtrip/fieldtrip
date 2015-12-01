function browse_audiovideo(cfg, data)

% BROWSE_AUDIOVIDEO reads and vizualizes the audio and/or video data
% corresponding to the EEG/MEG data segment that is passed into this
% function.

% this is a simple wrapper around the high-level function
cfg.interactive = 'no';
ft_audiovideobrowser(cfg, data);
