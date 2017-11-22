function test_bug1887

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_checkdata ft_datatype_raw ft_datatype_comp ft_datatype_timelock ft_componentanalysis ft_connectivityanalysis

% this contains raw data, 32 channels, 10 trials with nans
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1887.mat'));

for i=1:10
  data.trial{i} = randn(size(data.trial{i}));
end

cfg = [];
cfg.method = 'pca';
cfg.numcomponent = 5;
comp = ft_componentanalysis(cfg, data);

% create combined structures
cfg = [];
timelock = ft_timelockanalysis(cfg, comp);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
freq = ft_freqanalysis(cfg, comp);

%% simple checks
tmp = ft_checkdata(data, 'datatype', 'raw');

tmp = ft_checkdata(comp, 'datatype', 'comp');
tmp = ft_checkdata(comp, 'datatype', 'raw');
tmp = ft_checkdata(comp, 'datatype', 'raw+comp');

tmp = ft_checkdata(timelock, 'datatype', 'comp');
tmp = ft_checkdata(timelock, 'datatype', 'timelock');
tmp = ft_checkdata(timelock, 'datatype', 'timelock+comp');

tmp = ft_checkdata(freq, 'datatype', 'comp');
tmp = ft_checkdata(freq, 'datatype', 'freq');
tmp = ft_checkdata(freq, 'datatype', 'freq+comp');

%% this was problematic at a certain point
tmp = [];
tmp.topo = randn(3,3);
tmp.unmixing = randn(3,3);
tmp.topolabel = {'1', '2', '3'};
tmp = ft_checkdata(tmp, 'datatype', 'comp'); assert(~isempty(tmp));



%% the output should contain the component information
tmp = ft_checkdata(comp,     'datatype', 'raw+comp');      assert(isfield(tmp, 'topo'));
tmp = ft_checkdata(timelock, 'datatype', 'timelock+comp'); assert(isfield(tmp, 'topo'));
tmp = ft_checkdata(freq,     'datatype', 'freq+comp');     assert(isfield(tmp, 'topo'));

%% the output should NOT contain the component information
tmp = ft_checkdata(comp,     'datatype', 'raw');      assert(~isfield(tmp, 'topo'));
tmp = ft_checkdata(timelock, 'datatype', 'timelock'); assert(~isfield(tmp, 'topo'));
tmp = ft_checkdata(freq,     'datatype', 'freq');     assert(~isfield(tmp, 'topo'));

%% the output should NOT contain the "timeseries" information
tmp = ft_checkdata(comp,     'datatype', 'comp'); assert(~isfield(tmp, 'trial'));
tmp = ft_checkdata(timelock, 'datatype', 'comp'); assert(~isfield(tmp, 'avg'));
tmp = ft_checkdata(freq,     'datatype', 'comp'); assert(~isfield(tmp, 'powspctrm'));


%% the 2011 representation only supports raw timeseries data
tmp = ft_datatype_comp(comp,     'version', '2011'); assert(iscell(tmp.trial) && iscell(tmp.time));
tmp = ft_datatype_comp(timelock, 'version', '2011'); assert(iscell(tmp.trial) && iscell(tmp.time));

% the following does not work, as it requires time-frequency data for the conversion to raw
% tmp = ft_datatype_comp(freq, 'version', '2011'); assert(iscell(tmp.trial) && iscell(tmp.time)); 


%%
cfg = [];
cfg.layout = 'elec1020';
layout = ft_prepare_layout(cfg);

data = [];
data.label = layout.label;
data.trial = {};
for i=1:10
  data.time{i} = 1:1000;
  data.trial{i} = randn(length(data.label), 1000);
end

cfg = [];
cfg.method = 'pca';
cfg.numcomponent = 5;
comp = ft_componentanalysis(cfg, data);

% create combined structures
cfg = [];
timelock = ft_timelockanalysis(cfg, comp);

cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
freq = ft_freqanalysis(cfg, comp);

cfg = [];
cfg.layout = layout;
cfg.component = 1:3;
figure; ft_topoplotIC(cfg, comp);
figure; ft_topoplotIC(cfg, timelock);
figure; ft_topoplotIC(cfg, freq);

