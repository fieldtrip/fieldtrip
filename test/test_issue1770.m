function test_issue1770

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_preproc_dftfilter

% this is an issue related to numerical precision which is too strictly
% checked in ft_preproc_dftfilter (with spectral interpolation)

fs = 500;
data = [];
data.time{1} = -1:1/fs:2.4980; %-> observation: using linspace instead works fine

lnoise(1,:) = hanning(1750)'.*sin((2.*pi.*data.time{1}).*50);
lnoise(2,:) = hanning(1750)'.*sin((2.*pi.*data.time{1}).*100);

data.trial{1} = randn(2,1750)+lnoise;
data.label = {'a';'b'};
data.fsample = fs;

lineFreq = 50;

cfg              = [];
cfg.dftfilter    = 'yes'; %apply line noise filter with spectrum interpolation
cfg.dftfreq      = [lineFreq lineFreq*2]; %line noise and harmonic
cfg.dftreplace   = 'neighbour'; %spectral interpolation
cfg.dftbandwidth = [1 2]; %width of window to be interpolated
cfg.dftneighbourwidth = [2 2]; %width of window from which to interpolate
datafilt1 = ft_preprocessing(cfg, data);

cfg           = [];
cfg.dftfilter = 'yes';
cfg.dftfreq   = [lineFreq lineFreq*2];
datafilt2     = ft_preprocessing(cfg, data);

figure

subplot(2,2,1);plot(datafilt1.time{1}, data.trial{1}-datafilt1.trial{1}); ylim([-1.1 1.1]);xlabel('estimated linenoise method 1');
subplot(2,2,2);plot(datafilt2.time{1}, data.trial{1}-datafilt2.trial{1}); ylim([-1.1 1.1]);xlabel('estimated linenoise method 2');
subplot(2,2,3);plot(data.time{1}, lnoise); ylim([-1.1 1.1]);xlabel('simulated linenoise');


