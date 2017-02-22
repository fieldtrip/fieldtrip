% function inspect_ft_sourcemovie

% MEM 24gb
% WALLTIME 00:10:00

% TEST inspect_ft_sourcemovie
% TEST ft_sourcemovie ft_sourceanalysis ft_sourceinterpolate ft_volumenormalize ft_prepare_singleshell ft_prepare_leadfield qsubcellfun qsubfeval qsubget

% the frequency and source analysis is based on the tutorials

% use FieldTrip defaults instead of personal defaults
global ft_default;
ft_default = [];
ft_default.feedback = 'no';

% qsub is necessary, add fieldtrip/qsub to path
[v, p] = ft_version;
addpath(fullfile(p, 'qsub'));

%%

load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/timefrequencyanalysis/dataFIC.mat'))
load(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat'))

mri = ft_read_mri(dccnpath('/home/common/matlab/fieldtrip/data/Subject01.mri'));

cfg              = [];
cfg.output       = 'powandcsd';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:4:26;                         % changed from original
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -0.5:0.2:1.5;                   % changed from original
freqFIC          = ft_freqanalysis(cfg, dataFIC);

cfg = [];
vol = ft_prepare_singleshell(cfg, segmentedmri);

cfg                 = [];
cfg.grad            = freqFIC.grad;
cfg.vol             = vol;
cfg.reducerank      = 2;
cfg.normalize       = 'yes';
cfg.channel         = {'MEG','-MLP31', '-MLO12'};
cfg.grid.resolution = 1; % use a 3-D grid with a 1 cm resolution
cfg.grid.unit       = 'cm';
grid                = ft_prepare_leadfield(cfg);

%% do the source reconstruction

nfreq = length(freqFIC.freq);
ntime = length(freqFIC.time);

cfg = cell(nfreq, ntime);
for fbin=1:nfreq
  for tbin=1:ntime
    cfg{fbin,tbin}              = [];
    cfg{fbin,tbin}.latency      = freqFIC.time(tbin);
    cfg{fbin,tbin}.frequency    = freqFIC.freq(fbin);
    cfg{fbin,tbin}.method       = 'dics';
    cfg{fbin,tbin}.projectnoise = 'no';
    cfg{fbin,tbin}.grid         = grid;
    cfg{fbin,tbin}.vol          = vol;
    cfg{fbin,tbin}.lambda       = 0;
  end
end

memreq = 600*1024^2;  % requires ~400 MB
timreq = 45;          % requires ~15 seconds
source = qsubcellfun(@ft_sourceanalysis, cfg, repmat({freqFIC},nfreq,ntime), 'memreq', memreq, 'timreq', timreq);
clear cfg             % this one is very large

%% interpolate the coarse functional data onto the individuals anatomy
cfg             = [];
cfg.downsample  = 4;
cfg.parameter   = 'avg.pow';
memreq          = 500*1024^2;  % requires ~250 MB
timreq          = 24;          % requires ~8 seconds
sourceInt       = qsubcellfun(@ft_sourceinterpolate, repmat({cfg},nfreq,ntime), source, repmat({mri},nfreq,ntime), 'memreq', memreq, 'timreq', timreq);

if false
  cfg = [];
  cfg.funparameter = 'avg.pow';
  ft_sourceplot(cfg, sourceInt{1,1});
end

%% spatially normalize the individuals brain to the template MNI brain
cfg             = [];
cfg.coordsys    = 'ctf';
cfg.nonlinear   = 'no';
memreq          = 2*1024^3;    % requires ~1.3 GB
timreq          = 30;          % requires ~15 seconds
sourceIntNorm = qsubcellfun(@ft_volumenormalise, repmat({cfg}, nfreq, ntime), sourceInt, 'memreq', memreq, 'timreq', timreq);

if false
  cfg = [];
  cfg.funparameter = 'avg.pow';
  ft_sourceplot(cfg, sourceIntNorm{1,1});
end

%% interpolate the spatially normalized functional data on the cortical sheet
cortex = ft_read_headshape('cortex_8196.surf.gii');

cfg = [];
cfg.parameter = 'pow';
% memreq          = 150*1024^2;    % requires ~150 MB
% timreq          = 15;            % requires ~3 seconds
% sourceIntNormSurf = qsubcellfun(@ft_sourceinterpolate, repmat({cfg}, nfreq, ntime), sourceIntNorm, repmat({cortex}, nfreq, ntime), 'memreq', memreq, 'timreq', timreq);
sourceTFR
% for this one it is faster to run it locally than to go through the file system
sourceIntNormSurf = cellfun(@ft_sourceinterpolate, repmat({cfg}, nfreq, ntime), sourceIntNorm, repmat({cortex}, nfreq, ntime),'uniformoutput', false);

if false
  cfg = [];
  cfg.funparameter = 'pow';
  cfg.method = 'surface';
  ft_sourceplot(cfg, sourceIntNormSurf{1,1});
endsourceTFR

%%

npos = size(sourceIntNormSurf{1,1}.pos,1);

% construct a single source structure with one frequency and all time points
sourceER        = [];
sourceER.time   = freqFIC.time;
sourceER.pos    = sourceIntNormSurf{1,1}.pos;
sourceER.tri    = sourceIntNormSurf{1,1}.tri;
sourceER.dimord = 'pos_time';
sourceER.pow    = zeros(npos, ntime);
for tbin=1:ntime
  sourceER.pow(:,tbin) = sourceIntNormSurf{6,tbin}.pow;  % fbin 6 corresponds to 22 Hz
end

% construct a single source structure with all frequency and time points
sourceTFR        = [];
sourceTFR.time   = freqFIC.time;
sourceTFR.freq   = freqFIC.freq;
sourceTFR.pos    = sourceIntNormSurf{1,1}.pos;
sourceTFR.tri    = sourceIntNormSurf{1,1}.tri;
sourceTFR.dimord = 'pos_freq_time';
sourceTFR.pow    = zeros(npos, nfreq, ntime);
for fbin=1:nfreqsourceTFR
  for tbin=1:ntime
    sourceTFR.pow(:,fbin,tbin) = sourceIntNormSurf{fbin,tbin}.pow;
  end
end

%%

cfg = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg, sourceER);

%%

cfg = [];
cfg.funparameter = 'pow';
ft_sourcemovie(cfg, sourceTFR);
