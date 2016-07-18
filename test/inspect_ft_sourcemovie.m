function inspect_ft_sourcemovie

% MEM 1500mb
% WALLTIME 00:10:00

% TEST inspect_ft_sourcemovie
% TEST ft_sourcemovie ft_sourceanalysis ft_prepare_singleshell ft_prepare_leadfield
% TEST qsubcellfun qsubfeval qsubget

% the frequency and source analysis is based on the tutorials

% NOT FINISHED YET, SEE BELOW
return

global ft_default;
ft_default.feedback = 'no';

% qsub is necessary, add fieldtrip/qsub to path
addpath('/home/common/matlab/fieldtrip/qsub/')

load /home/common/matlab/fieldtrip/data/ftp/tutorial/timefrequencyanalysis/dataFIC.mat
load /home/common/matlab/fieldtrip/data/ftp/tutorial/beamformer/segmentedmri.mat

mri = ft_read_mri(fullfile(homedir, 'common', 'matlab', 'fieldtrip', 'data', 'Subject01.mri'));

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
cfg.channel         = {'MEG','-MLP31', '-MLO12'};
cfg.grid.resolution = 1; % use a 3-D grid with a 1 cm resolution
[grid]              = ft_prepare_leadfield(cfg);

nfreq = length(freqFIC.freq);
ntime = length(freqFIC.time);

clear cfg
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

% interpolate the coarse functional data onto the individuals anatomy
cfg             = [];
cfg.downsample  = 4;
cfg.parameter   = 'avg.pow';
memreq          = 500*1024^2;  % requires ~250 MB
timreq          = 24;          % requires ~8 seconds
sourceInt       = qsubcellfun(@ft_sourceinterpolate, repmat({cfg},nfreq,ntime), source, repmat({mri},nfreq,ntime), 'memreq', memreq, 'timreq', timreq);

% spatially normalize the individuals brain to the template brain
cfg             = [];
cfg.coordsys    = 'ctf';
cfg.nonlinear   = 'no';
cfg.downsample  = 4;
memreq          = 60*1024^2;   % requires ~30 MB
timreq          = 18;          % requires ~6 seconds
sourceIntNorm = qsubcellfun(@ft_volumenormalise, repmat({cfg}, nfreq, ntime), sourceInt, 'memreq', memreq, 'timreq', timreq);

% interpolate the functional data on the cortical sheet
% FIXME to be continued
error('right now I do not know how to do the interpolation on the cortical sheet')

npos = size(sourceIntNorm{1,1}.pos,1);

% construct a single source structure with one frequency and all time points
sourceER        = [];
sourceER.time   = freqFIC.time;
sourceER.pos    = sourceIntNorm{1,1}.pos;
sourceER.dim    = sourceIntNorm{1,1}.dim;
sourceER.dimord = 'pos_time'
sourceER.pow    = zeros(npos, ntime);
for tbin=1:ntime
sourceER.pow(:,tbin) = sourceIntNorm{6,tbin}.avg.pow;  % fbin 6 corresponds to 22 Hz
end

% construct a single source structure with all frequency and time points
sourceTFR        = [];
sourceTFR.time   = freqFIC.time;
sourceTFR.freq   = freqFIC.freq;
sourceTFR.pos    = sourceIntNorm{1,1}.pos;
sourceTFR.dim    = sourceIntNorm{1,1}.dim;
sourceTFR.dimord = 'pos_freq_time'
sourceTFR.pow    = zeros(npos, nfreq, ntime);
for fbin=1:nfreq
for tbin=1:ntime
sourceTFR.pow(:,fbin,tbin) = sourceIntNorm{fbin,tbin}.avg.pow;
end
end

cfg = [];
ft_sourcemovie(cfg, sourceER);

cfg = [];
ft_sourcemovie(cfg, sourceTFR);

