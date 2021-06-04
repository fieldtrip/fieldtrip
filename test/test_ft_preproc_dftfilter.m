function test_ft_preproc_dftfilter

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some other instances, testing ft_preproc_dftfilter directly
tim = (0:999)./1000;
dat = randn(1, 1000) + hanning(1000)'.*sin(2.*pi.*(tim).*50);
filt = ft_preproc_dftfilter(dat, 1000, 50, 'dftreplace', 'neighbour');
figure;plot(tim, dat-filt); hold on;plot(tim, hanning(1000)'.*sin(2.*pi.*(tim).*50));

tim = (0:1000)./1000;
dat = randn(1, 1001) + hanning(1001)'.*sin(2.*pi.*(tim).*50);
filt = ft_preproc_dftfilter(dat, 1000, 50, 'dftreplace', 'neighbour');
figure;plot(tim, dat-filt); hold on;plot(tim, hanning(1001)'.*sin(2.*pi.*(tim).*50));


dat = randn(3, 1001) + [hanning(1001)'.*sin(2.*pi.*(tim).*53) ; hanning(1001)'.*sin(2.*pi.*(tim).*79 - 0.025) ; hanning(1001)'.*sin(2.*pi.*(tim).*127 + 0.002)];
filt = ft_preproc_dftfilter(dat, 1000, [53, 79, 127], 'dftreplace', 'neighbour', 'dftneighbourwidth', [1 1 1]);
figure; 
subplot(2,2,1); plot(tim, dat(1,:)-filt(1,:)); hold on;plot(tim, hanning(1001)'.*sin(2.*pi.*(tim).*53));
subplot(2,2,2); plot(tim, dat(2,:)-filt(2,:)); hold on;plot(tim, hanning(1001)'.*sin(2.*pi.*(tim).*79 - 0.025));
subplot(2,2,3); plot(tim, dat(3,:)-filt(3,:)); hold on;plot(tim, hanning(1001)'.*sin(2.*pi.*(tim).*127 + 0.002));

tim = (0:1000)./678.253;
krn = hanning(1000)';
dat = randn(1, 1001) + ([krn(1:500) ones(1,501)]).*sin(2.*pi.*(tim).*50);
filt = ft_preproc_dftfilter(dat, 678.253, 50, 'dftreplace', 'neighbour', 'dftneighbourwidth', 4, 'dftbandwidth', 2);
figure;plot(tim, dat-filt); hold on;plot(tim, ([krn(1:500) ones(1,501)]).*sin(2.*pi.*(tim).*50));

