function test_old_connectivityanalysis

% MEM 1gb
% WALLTIME 00:10:00


% this script tests the functionality of connectivityanalysis

is_octave=~ft_platform_supports('matlabversion',1,inf);
if is_octave
  % TODO: figure out what the issue is here
  % Note: this seems hard to reproduce on non-Travis system
  reason=sprintf('%s can crash Octave on travis and has been distabled',...
                    mfilename());
  moxunit_throw_test_skipped_exception(reason)
end

% first create some data
cfg             = [];
cfg.method      = 'linear_mix'; 
cfg.ntrials     = 200;
cfg.nsignal     = 2;
cfg.triallength = 1;
cfg.fsample     = 200;
%cfg.fsample     = 1000;
%cfg.bpfilter    = 'yes';
%cfg.bpfreq      = [15 25];
%cfg.mix         = [1 0 1; 0 1 1];
%cfg.delay       = [0 0 5; 0 0 0];
cfg.method      = 'ar';
cfg.params(:,:,1) = [ 0.8 0.5; 
                      0   0.9];
cfg.params(:,:,2) = [-0.5    0; 
                        0 -0.8]; 
cfg.noisecov      = [0.3 0; 
                       0 1];
data            = ft_connectivitysimulation(cfg);

% do mvaranalysis
cfgm       = [];
cfgm.order = 5;
cfgm.toolbox = 'bsmart';
mdata      = ft_mvaranalysis(cfgm, data);
cfgfm      = [];
cfgfm.method = 'mvar';
mfreq      = ft_freqanalysis(cfgfm, mdata);

% freqanalysis
cfgf           = [];
cfgf.method    = 'mtmfft';
cfgf.output    = 'fourier';
cfgf.tapsmofrq = 2;
freq           = ft_freqanalysis(cfgf, data);
fd             = ft_freqdescriptives([], freq);
%freqx          = freq2transfer([], freq);

%cfgsf.channelcmb = {'all' 'all'};
%freq2x2          = freq2transfer(cfgsf, freq);

% connectivityanalysis
cfgc           = [];
cfgc.method    = 'coh';
c1             = ft_connectivityanalysis(cfgc, freq);
c1m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'plv';
c2             = ft_connectivityanalysis(cfgc, freq);
c2m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'csd';
cfgc.complex   = 'angle';
c3             = ft_connectivityanalysis(cfgc, freq);
c3m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'psi';
cfgc.complex   = 'abs';
cfgc.bandwidth = 4;
c4             = ft_connectivityanalysis(cfgc, freq);
c4m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'granger';
c5             = ft_connectivityanalysis(cfgc, freq);
c5m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.granger.sfmethod = 'bivariate';
c5b            = ft_connectivityanalysis(cfgc, freq);
cfgc = rmfield(cfgc, 'granger');

cfgc.method    = 'pdc';
c6             = ft_connectivityanalysis(cfgc, freq);
c6m            = ft_connectivityanalysis(cfgc, mfreq);

cfgc.method    = 'dtf';
c7             = ft_connectivityanalysis(cfgc, freq);
c7m            = ft_connectivityanalysis(cfgc, mfreq);

cfgc.method    = 'instantaneous_causality';
c8             = ft_connectivityanalysis(cfgc, freq);
c8m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.granger.sfmethod = 'bivariate';
c8c            = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.sfmethod = 'multivariate';


cfgc.method    = 'total_interdependence';
c9             = ft_connectivityanalysis(cfgc, freq);
c9m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.granger.sfmethod = 'bivariate';
c9c            = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.sfmethod = 'multivariate';


%--------------------------------------------------------
% now make 3 channels with no direct link between 1 and 2
% create some data
cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
%cfg.method      = 'linear_mix'; 
%cfg.bpfilter    = 'yes';
%cfg.bpfreq      = [15 25];
%cfg.mix         = [1 0 1; 0 1 1;0 0 2];
%cfg.delay       = [0 0 5; 0 0 0;0 0 0];
cfg.method      = 'ar';
cfg.params(:,:,1) = [ 0.8 0   0; 
                      0   0.9 0.5;
                      0.4 0   0.5];
cfg.params(:,:,2) = [-0.5    0  0; 
                        0 -0.8  0; 
                        0    0 -0.2];
cfg.noisecov      = [0.3 0 0;
                       0 1 0;
                       0 0 0.2];

data            = ft_connectivitysimulation(cfg);

% do mvaranalysis
cfgm       = [];
cfgm.order = 5;
cfgm.toolbox = 'bsmart';
mdata      = ft_mvaranalysis(cfgm, data);
cfgfm      = [];
cfgfm.method = 'mvar';
mfreq      = ft_freqanalysis(cfgfm, mdata);

% freqanalysis
cfgf           = [];
cfgf.method    = 'mtmfft';
cfgf.output    = 'fourier';
cfgf.tapsmofrq = 2;
freq           = ft_freqanalysis(cfgf, data);
fd             = ft_freqdescriptives([], freq);

cfgfs            = [];
cfgfs.channelcmb = {'all' 'all'};

% connectivityanalysis
cfgc           = [];
cfgc.method    = 'coh';
c1             = ft_connectivityanalysis(cfgc, freq);
c1m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'plv';
c2             = ft_connectivityanalysis(cfgc, freq);
c2m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'csd';
cfgc.complex   = 'angle';
c3             = ft_connectivityanalysis(cfgc, freq);
c3m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'psi';
cfgc.complex   = 'abs';
cfgc.bandwidth = 4;
c4             = ft_connectivityanalysis(cfgc, freq);
c4m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'granger';
c5             = ft_connectivityanalysis(cfgc, freq);
c5m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.granger.sfmethod = 'bivariate';
c5b            = ft_connectivityanalysis(cfgc, freq);
cfgc = rmfield(cfgc, 'granger');

cfgc.method    = 'pdc';
c6             = ft_connectivityanalysis(cfgc, freq);
c6m            = ft_connectivityanalysis(cfgc, mfreq);

cfgc.method    = 'dtf';
c7             = ft_connectivityanalysis(cfgc, freq);
c7m            = ft_connectivityanalysis(cfgc, mfreq);

cfgc.method    = 'instantaneous_causality';
c8             = ft_connectivityanalysis(cfgc, freq);
c8m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'total_interdependence';
c9             = ft_connectivityanalysis(cfgc, freq);
c9m            = ft_connectivityanalysis(cfgc, mfreq);



cfgc             = [];
cfgc.partchannel = 'signal003';
cfgc.method      = 'csd';
cfgc.complex     = 'complex';
freqp            = ft_connectivityanalysis(cfgc, freq);
cfgc.method      = 'granger';
c5p              = ft_connectivityanalysis(cfgc, freqp);
cfgc.method      = 'instantaneous_causality';
c8p              = ft_connectivityanalysis(cfgc, freqp);



