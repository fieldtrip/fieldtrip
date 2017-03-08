function test_bug2664

% MEM 2gb
% WALLTIME 00:10:00

% TEST ft_sourceanalysis ft_checkdata

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];

%% create some data
data = [];
data.label = {'a' 'b' 'c'};
data.fsample = 256;
data.trial = {randn(3,1024)};
data.time = {1/256:1/256:4};

% not sure what this is for but the original bug report had this field
data.interpolatedElectrodes = {randn(3,1024)};

%% attempt source analysis MNE

cfg = struct;
cfg.method = 'mne';

% the following will lead to an error because the fields are empty, but the
% error we are testing for is a different one ("Reference to non-existent
% field 'topo'." )
cfg.elec = [];
cfg.grid = [];
cfg.vol = [];
cfg.hdmfile = [];

cfg.rawtrial = 'yes';
cfg.mne.lambda = '5%';
cfg.keepfilter = 'yes';
cfg.rawtrial = 'no';
cfg.singletrial = 'no';
cfg.keeptrials = 'yes';

try
  source = ft_sourceanalysis(cfg, data);
catch err
  if strcmp(err.message, 'Reference to non-existent field ''topo''.')
    rethrow(err);
  end % ignore all other errors for this test
end

end
