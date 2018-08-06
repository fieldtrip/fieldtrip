function test_bug1571

% MEM 3gb
% WALLTIME 00:10:00

% TEST ft_preprocessing ft_preproc_dftfilter dftfilter ft_freqanalysis ft_singleplotER

% Philipp Hintze wrote:
% The issue with the continuous data is the following: I tried using
% [this script]
% but looking at averaged timelocked data for the subject, the data looks
% identical to the result I get using no filter at all, i.e., contaminated by
% strong line noise. Including padding makes no difference either, therefore
% the brackets.

dataset = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1571/20101230_1010.cnt');

cfg                     = [];
cfg.channel      = {'EEG'};
cfg.datafile     = dataset;
cfg.headerfile   = dataset;
cfg.dataset      = dataset;
cfg.continuous   = 'yes';

unfilteredContinuousData  = ft_preprocessing(cfg);

cfg.dftfilter           = 'yes';
filteredContinuousData  = ft_preprocessing(cfg);

assert(~isequal(filteredContinuousData.trial{1}(:), unfilteredContinuousData.trial{1}(:)))

cfg = [];
cfg.channel = {'Fp1'};
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foi = 1:100;
fftFilt = ft_freqanalysis(cfg, filteredContinuousData);
figure;
subplot(2,1,1)
ft_singleplotER(cfg, fftFilt)

cfg = [];
cfg.channel = {'Fp1'};
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foi = 1:100;
fftUnfilt = ft_freqanalysis(cfg, unfilteredContinuousData);
subplot(2,1,2)
ft_singleplotER(cfg, fftUnfilt)

%% synthetic signal with wrong sampling rate

fsample = 512;

data = {};
data.fsample = fsample;
data.label = {'1', '2', '3'};
nchans = numel(data.label);
real_fsample = 512.01;
data.time{1} = 0:1/real_fsample:100;
data.trial{1} = rand(nchans, numel(data.time{1}));

cfg                     = [];
cfg.continuous          = 'yes';
cfg.dftfilter           = 'yes';
filteredData  = ft_preprocessing(cfg, data);

assert(~isequal(data.trial{1}, filteredData.trial{1}(:)))

fprintf('done.\n');


