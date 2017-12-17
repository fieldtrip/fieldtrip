function test_firwsfiltering

% MEM 2048mb
% WALLTIME 00:10:00

% TEST ft_preprocessing filter_with_correction ft_preproc_bandpassfilter
% TEST ft_preproc_highpassfilter ft_preproc_lowpassfilter firws

%% create some data

data = [];
data.label = {'a','b','c'};

nTri = 50;
data.time = cell(nTri,1);
data.trial = cell(nTri,1);
for k = 1:nTri
  data.time{k} = 0.001:0.001:1;
  data.trial{k} = randn(3,1000);
end

% introduce a time axis different from the rest
data.time{10} = 0.001:0.001:3;
data.trial{10} = randn(3,3000);

data.time{1} = 0.001:0.001:8;
data.trial{1} = randn(3,8000);

%% filter it

cfg = [];
cfg.lpfilter = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq = 20;

% plot the filter response
cfg.plotfiltresp = 'yes';

datfilt = ft_preprocessing(cfg, data);

%%

cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfilttype = 'firws';
cfg.bpfreq = [20 30];

% plot the filter response
cfg.plotfiltresp = 'yes';

datfilt = ft_preprocessing(cfg, data);

%%

cfg = [];
cfg.hpfilter = 'yes';
cfg.hpfilttype = 'firws';
cfg.hpfreq = 20;

% plot the filter response
cfg.plotfiltresp = 'yes';

datfilt = ft_preprocessing(cfg, data);

%%

cfg = [];
cfg.bsfilter = 'yes';
cfg.bsfilttype = 'firws';
cfg.bsfreq = [20 30];

% plot the filter response
cfg.plotfiltresp = 'yes';

datfilt = ft_preprocessing(cfg, data);

end
