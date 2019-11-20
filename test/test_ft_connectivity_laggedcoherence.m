function test_ft_connectivity_laggedcoherence

% MEM 4gb
% WALLTIME 00:10:00
% DEPENDENCY ft_connectivityanalysis ft_laggedcoherence fourierspctrm2lcrsspctrm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2951
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --------------------------------------------------------
% first create some data, this is the same as in test_ft_connectivityanalysis
% make 3 channels with no direct link between 1 and 2
cfgs             = [];
cfgs.ntrials     = 500;
cfgs.triallength = 1;
cfgs.fsample     = 200;
cfgs.nsignal     = 3;
cfgs.method      = 'ar';

cfgs.params(:,:,1) = [ 0.8  0.0  0.0;
                       0.0  0.9  0.5;
                       0.4  0.0  0.5];

cfgs.params(:,:,2) = [-0.5  0.0  0.0;
                       0.0 -0.8  0.0;
                       0.0  0.0 -0.2];

cfgs.noisecov      = [ 0.3  0.0  0.0;
                       0.0  1.0  0.0;
                       0.0  0.0  0.2];

data = ft_connectivitysimulation(cfgs);

%% --------------------------------------------------------
% freqanalysis
cfgf           = [];
cfgf.method    = 'mtmconvol';
cfgf.output    = 'fourier';
cfgf.taper     = 'hanning';
cfgf.toi       = 0:0.01:1;
cfgf.foi       = 2:2:70;
cfgf.tapsmofrq = 2   * ones(size(cfgf.foi));
cfgf.t_ftimwin = 0.5 * ones(size(cfgf.foi));
freq           = ft_freqanalysis(cfgf, data);

%% --------------------------------------------------------
% connectivityanalysis
cfgc           = [];
cfgc.method    = 'laggedcoherence';
c1             = ft_connectivityanalysis(cfgc, freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% see https://github.com/fieldtrip/fieldtrip/issues/1217 and https://github.com/fieldtrip/fieldtrip/pull/1233
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1217.mat'));

% do a spectral decomposition first
cfg = [];
cfg.method = 'mtmconvol';
cfg.output = 'fourier';
cfg.foi    = 2:2:100;
cfg.t_ftimwin = ones(1,numel(cfg.foi))./2;
cfg.taper  = 'hanning';
cfg.toi    = -0.5:0.05:1.5;
freq = ft_freqanalysis(cfg, data);

cfg = [];
cfg.method = 'laggedcoherence';
cfg.laggedcoherence.lags = 0.5;
lcoh = ft_connectivityanalysis(cfg, freq);

% try it on a single frequency
cfg = [];
cfg.frequency = 10;
freq2 = ft_selectdata(cfg, freq);

cfg = [];
cfg.method = 'laggedcoherence';
cfg.laggedcoherence.lags = 0.5;
lcoh2 = ft_connectivityanalysis(cfg, freq2);
