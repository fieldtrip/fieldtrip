function test_ft_preproc_medianfilter

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_preproc_medianfilter

% tic
inspect_results = false; %true if visual inspection required

nconfig = 2; %number of different configurations to test

fs = 500;
nchan = 32;
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
cfg.medianfilter    = 'yes'; %apply line noise filter with spectrum interpolation
datafilt{1} = ft_preprocessing(cfg, data);

choices = num2cell(2:6);
for i = 1:length(choices)
    cfg.medianfiltord = choices{i};
    datafilt{2} = ft_preprocessing(cfg, data);
end

if inspect_results
    figure;
    for i = 1:nconfig
        subplot(ceil(nconfig/2),2,i);plot(datafilt{i}.time{1}, datafilt{i}.trial{1}(1,:)); ylim([-2.1 2.1]);xlabel(['config ' num2str(i)]);
    end
end

clear datafilt
% toc