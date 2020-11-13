function test_ricoh

% MEM 6gb
% WALLTIME 01:00:00
% DEPENDENCY hasricoh isricohmegfile read_ricoh_data read_ricoh_event read_ricoh_header ricoh2grad

%% Test function for Ricoh *.con (Third-party exported) data
%
% Recording data of a session of AEF for right-ear stimulation composed of 100 trials
% This is an "exported .con" data that includes full of recording information:
% MEG data, EEG data, digitized points, acquisition conditions, etc.

% filename = 'rightearAEF_export.con';
filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/ricoh160/rightearAEF_export.con');

if ~ft_filetype(filename, 'ricoh_con')
  error('definition error of ricoh_con');
end

%% ft_read_header
hdr = ft_read_header(filename)

if length(ft_channelselection('MEG', hdr.label))~=160
  error('did not select all MEG channels');
end

if length(ft_channelselection('eeg', hdr.label))~=42
  error('did not select all EEG channels');
end


if ~ft_senstype(hdr, 'yokogawa160')
  error('hdr ~= yokogawa160');
end

%% ft_preprocessing
cfg = [];
cfg.dataset  = filename;
data = ft_preprocessing(cfg)

%% ft_read_event
event = ft_read_event(filename, 'threshold', 1.6)

%% Define trials
cfg = [];
cfg.dataset          = filename;
cfg.trialfun  = 'ft_trialfun_general';
cfg.trialdef.eventtype  = 'triginfo';
cfg.trialdef.eventvalue = 'AEF';
cfg.trialdef.prestim    = 0.1; %sec
cfg.trialdef.poststim   = 0.5; %sec
cfg = ft_definetrial(cfg);
trl = cfg.trl;

if length(trl(:,1))~=100
  error('could not define trials appropriately');
end

%% ft_read_headshape
headshape = ft_read_headshape(filename)

%% ft_read_sens
grad = ft_read_sens(filename)

if length(ft_chanunit(grad,'T'))~=160
  error('did not select all MEG channels');
end
