% this script tests the functionality of FT_CONNECTIVITYANALYSIS
% on frequency domain channel data

% apart from using FT_CONNECTIVITYANALYSIS
% it also relies on FT_CONNECTIVITYSIMULATION, FT_FREQANALYSIS, FT_MVARANALYSIS

% first create some data
%--------------------------------------------------------
% make 3 channels with no direct link between 1 and 2
cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
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
cfgc.sfmethod  = 'bivariate';
c5b            = ft_connectivityanalysis(cfgc, freq);
cfgc           = rmfield(cfgc, 'sfmethod');
cfgc.method    = 'pdc';
c6             = ft_connectivityanalysis(cfgc, freq);
c6m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'dtf';
c7             = ft_connectivityanalysis(cfgc, freq);
c7m            = ft_connectivityanalysis(cfgc, mfreq);
%cfgc.method    = 'instantaneous_causality';
%c8             = ft_connectivityanalysis(cfgc, freq);
%c8m            = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'total_interdependence';
c9             = ft_connectivityanalysis(cfgc, freq);
c9m            = ft_connectivityanalysis(cfgc, mfreq);



cfgc             = [];
cfgc.partchannel = 'signal003';
cfgc.method      = 'csd';
cfgc.complex     = 'complex';
freqp            = ft_connectivityanalysis(cfgc, freqx);
freqxp           = freq2transfer([], freqp);
cfgc.method      = 'granger';
c5p              = ft_connectivityanalysis(cfgc, freqxp);
cfgc.method      = 'instantaneous_causality';
c8p              = ft_connectivityanalysis(cfgc, freqxp);



