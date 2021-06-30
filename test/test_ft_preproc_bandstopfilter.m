function test_ft_preproc_bandstopfilter

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_preproc_bandstopfilter

% tic
inspect_results = false; %true if visual inspection required

nconfig = 8; %number of different configurations to test

fs = 500;
nchan = 32;
stopband = [49,51 ; 98,102];
start_time = -1; %s
end_time = 2.5; %s
nsamples = (end_time - start_time) * fs;
data = [];
data.time{1} = linspace(start_time, end_time, nsamples);

data.trial{1} = randn(nchan,nsamples);
data.label = cellstr(num2str((1:nchan).'));
data.fsample = fs;

datafilt = cell(1,nconfig);

cfg              = [];
cfg.bsfilter    = 'yes'; %apply line noise filter with spectrum interpolation
cfg.bsfreq      = stopband;
datafilt{1} = ft_preprocessing(cfg, data);

choices = num2cell(2:6);
for i = 1:length(choices)
    cfg.bsfiltord = choices{i};
    datafilt{2} = ft_preprocessing(cfg, data);
end

choices = {'firws' , 'fir' , 'firls'};
for i = 1:length(choices)
    cfg.bsfilttype = choices{i};
    datafilt{3} = ft_preprocessing(cfg, data);
end

choices = {'onepass' , 'onepass-reverse' , 'onepass-zerophase'};
for i = 1:length(choices)
    cfg.bsfiltdir = choices{i};%= filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse' or 'onepass-zerophase' (default for firws) or 'onepass-minphase' (firws, non-linear!)
    datafilt{4} = ft_preprocessing(cfg, data);
end

choices = {'reduce', 'split'};
for i = 1:length(choices)
    cfg.bsinstabilityfix = choices{i};
    datafilt{5} = ft_preprocessing(cfg, data);
end

choices = num2cell(1:5);
for i = 1:length(choices)
    cfg.bsfiltdf = choices{i};
    datafilt{6} = ft_preprocessing(cfg, data);
end

choices = {'hann' , 'blackman' , 'kaiser'};
for i = 1:length(choices)
    cfg.bsfiltwintype = choices{i};
    datafilt{7} = ft_preprocessing(cfg, data);
end

choices = num2cell(logspace(-3,0,4));
for i = 1:length(choices)
    cfg.bsfiltdev = choices{i};
    datafilt{8} = ft_preprocessing(cfg, data);
end


if inspect_results
    figure;
    for i = 1:nconfig
        subplot(ceil(nconfig/2),2,i);plot(datafilt{i}.time{1}, datafilt{i}.trial{1}(1,:)); ylim([-2.1 2.1]);xlabel(['config ' num2str(i)]);
    end
end

clear datafilt
% toc