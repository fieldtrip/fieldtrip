function test_ft_preproc_baselinecorrect

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_preproc_baselinecorrect

% tic
inspect_results = false; %true if visual inspection required

nconfig = 4; %number of different configurations to test

fs = 500;
nchan = 32;
start_time = -1; %s
end_time = 2.5; %s
nsamples = (end_time - start_time) * fs;
timevect = linspace(start_time, end_time, nsamples)';
% label = cellstr(num2str((1:nchan).'));

begsample = 10; %time (in seconds) of the begin sample for the baseline estimate
endsample = 510; %time (in seconds) of the end sample for the baseline estimate
offset = 1.5;
data = randn(nchan,nsamples)+offset;

dataproc = cell(1,nconfig);
dataproc{1} = ft_preproc_baselinecorrect(data);
dataproc{2} = ft_preproc_baselinecorrect(data,begsample);
dataproc{3} = ft_preproc_baselinecorrect(data,begsample,endsample);

data2 = [];
data2.time{1} = timevect;
data2.trial{1} = data;
data2.label = cellstr(num2str((1:nchan).'));
data2.fsample = fs;

% test in ft_preprocessing
cfg = [];
cfg.demean = 'yes';
cfg.baselinewindow = [begsample,endsample];
% dataproc{4} = preproc(data, label, timevect, cfg);
dataproc{4} = ft_preprocessing(cfg,data2);

if inspect_results
    figure;
    for i = 1:nconfig-1
        subplot(ceil(nconfig/2),2,i);plot(timevect, dataproc{i}(1,:)-data(1,:)); ylim([-2.1 2.1]);xlabel(['config ' num2str(i)]);
    end
    i = 4;
    subplot(ceil(nconfig/2),2,i);plot(timevect, mean(dataproc{i}.trial{1})-mean(data)); ylim([-2.1 2.1]);xlabel(['config ' num2str(i)]);
end

clear dataproc
% toc

