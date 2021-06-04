function test_issue1770

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_preproc_dftfilter

% this is an issue related to numerical precision which is too strictly
% checked in ft_preproc_dftfilter (with spectral interpolation)

fs = 500;
data = [];
data.time{1} = -1:1/fs:2.4980; %-> observation: using linspace instead works fine
data.label = {'a';'b'};
data.fsample = fs;

lineFreq = 50;

cfg              = [];
cfg.dftfilter    = 'yes'; %apply line noise filter with spectrum interpolation
cfg.dftfreq      = [lineFreq lineFreq*2]; %line noise and harmonic
cfg.dftreplace   = 'neighbour'; %spectral interpolation
cfg.dftbandwidth = [1 2]; %width of window to be interpolated
cfg.dftneighbourwidth = [2 2]; %width of window from which to interpolate
datafilt = ft_preprocessing(cfg, data);
