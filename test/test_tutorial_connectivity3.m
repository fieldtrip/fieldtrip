function test_tutorial_connectivity3

% TEST test_tutorial_connectivity3
% TEST ft_timelockanalysis ft_sourceanalysis ft_connectivityanalysis

% This is the third section of the connectivity tutorial, which
% starts with the CMC dataset, extracts a virtual channel and performs
% connectivity analysis on the virtual channel time series.

load source
[maxval, maxindx] = max(source.avg.coh);
maxpos = source.pos(maxindx,:)
clear source

load data

cfg               = [];
cfg.covariance    = 'yes';
cfg.channel       = 'MEG';
% cfg.keeptrials    = 'no';
cfg.vartrllength  = 2;
timelock          = ft_timelockanalysis(cfg, data);

cfg             = [];
cfg.method      = 'lcmv';
cfg.hdmfile     = 'SubjectCMC.hdm';
cfg.inwardshift = 1;
cfg.grid.pos    = maxpos;
cfg.keepfilter  = 'yes';
source          = ft_sourceanalysis(cfg, timelock);

chansel = ft_channelselection('MEG', data.label); % find the names
chansel = match_str(data.label, chansel);         % find the indices

beamformer = source.avg.filter{1};

virtual = [];
virtual.label = {'x', 'y', 'z'};
virtual.time = data.time;
for i=1:length(data.trial)
  virtual.trial{i} = beamformer * data.trial{i}(chansel,:);
end

timeseries = cat(2, data.trial{i});

[u, s, v] = svd(timeseries, 'econ');


