function test_tutorial_connectivity3(datadir)

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_timelockanalysis ft_sourceanalysis ft_connectivityanalysis ft_prepare_sourcemodel headsurface

% This is the third section of the connectivity tutorial, which
% starts with the CMC dataset, extracts a virtual channel and performs
% connectivity analysis on the virtual channel time series.

global ft_default;
ft_default.feedback = 'no';
ft_default.checkconfig = 'loose';

if nargin==0
  % this is where the data should be located
  datadir = dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/connectivity');
end

load(fullfile(datadir, 'source.mat'));

[maxval, maxindx] = max(source.avg.coh);
maxpos = source.pos(maxindx,:);

load(fullfile(datadir,'data.mat'));

%% compute the beamformer filter
cfg                   = [];
cfg.covariance        = 'yes';
cfg.channel           = 'MEG';
cfg.vartrllength      = 2;
cfg.covariancewindow  = 'all';
timelock              = ft_timelockanalysis(cfg, data);

cfg             = [];
cfg.method      = 'lcmv';
cfg.hdmfile     = fullfile(datadir,'SubjectCMC.hdm');
cfg.grid.pos    = maxpos;
cfg.keepfilter  = 'yes';
source          = ft_sourceanalysis(cfg, timelock);

%% construct the 3-D virtual channel at the location of interest
beamformer = source.avg.filter{1};

chansel = ft_channelselection('MEG', data.label); % find the names
chansel = match_str(data.label, chansel);         % find the indices

sourcedata = [];
sourcedata.label = {'x', 'y', 'z'};
sourcedata.time = data.time;
for i=1:length(data.trial)
  sourcedata.trial{i} = beamformer * data.trial{i}(chansel,:);
end

cfg = [];
cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
ft_databrowser(cfg, sourcedata);

%% construct a single virtual channel in the maximum power orientation
timeseries = cat(2, sourcedata.trial{:});

[u, s, v] = svd(timeseries, 'econ');

% whos u s v
%   Name           Size              Bytes  Class     Attributes
% 
%   s              3x3                  72  double              
%   u              3x3                  72  double              
%   v         196800x3             4723200  double            
  
% this is equal to the first column of matrix V, apart from the scaling with s(1,1)
timeseriesmaxproj = u(:,1)' * timeseries;

virtualchanneldata = [];
virtualchanneldata.label = {'cortex'};
virtualchanneldata.time = data.time;
for i=1:length(data.trial)
  virtualchanneldata.trial{i} = u(:,1)' * beamformer * data.trial{i}(chansel,:);
end

%% combine the virtual channel with the two EMG channels
cfg = [];
cfg.channel = 'EMG';
emgdata = ft_selectdata(cfg, data);

cfg = [];
combineddata = ft_appenddata(cfg, virtualchanneldata, emgdata);

% save combineddata combineddata

%% compute the spectral decomposition
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 100];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channel    = {'cortex' 'EMGlft' 'EMGrgt'};
freq    = ft_freqanalysis(cfg, combineddata);

cfg = [];
cfg.method = 'coh';
coherence = ft_connectivityanalysis(cfg, freq);

cfg = [];
cfg.zlim = [0 0.2];
figure
ft_connectivityplot(cfg, coherence);
title('coherence')

figure
plot(coherence.freq, squeeze(coherence.cohspctrm(1,2,:)))
title(sprintf('connectivity between %s and %s', coherence.label{1}, coherence.label{2}));
xlabel('freq (Hz)')
ylabel('coherence')





