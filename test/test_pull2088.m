function test_pull2088(filename)

% WALLTIME 00:20:00
% MEM 3gb
% DEPENDENCY ft_definetrial ft_trialfun_general
% DATA private

if nargin<1
  filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeglab/sub-001-eeg-sub-001_task-P300_run-2_eeg.set');
else
  % use the filename from the input, if not on the filesystem at dccn
end
  
cfg              = [];
cfg.dataset      = filename;
cfg.trialdef.eventtype  = 'trigger';
cfg.trialdef.eventvalue = { 'standard' 'oddball_with_reponse' };
cfg.trialdef.prestim        = 1; % in seconds
cfg.trialdef.poststim       = 2; % in seconds
cfg = ft_definetrial(cfg);
cfg.trl = cfg.trl(1:50,:);

% Baseline-correction options
cfg.demean          = 'yes';
cfg.baselinewindow  = [-0.3 0];
data = ft_preprocessing(cfg);
assert(isfield(data, 'trialinfo'));
