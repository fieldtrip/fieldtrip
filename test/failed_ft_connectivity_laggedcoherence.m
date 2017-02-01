function failed_ft_connectivity_laggedcoherence

% WALLTIME 00:10:00
% MEM 3gb

% see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=2951

% TEST test_ft_connectivity_laggedcoherence
% TEST ft_connectivity_laggedcoherence ft_connectivityanalysis

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

warning('the following will give an error because the implementation is still incomplete');

c1             = ft_connectivityanalysis(cfgc, freq);

