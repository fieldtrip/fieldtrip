function test_ft_connectivityanalysis

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_connectivityanalysis ft_connectivity_granger ft_connectivity_corr ft_connectivity_psi ft_mvaranalysis ft_connectivitysimulation ft_freqanalysis ft_connectivity_pdc ft_connectivity_dtf ft_connectivity_csd2transfer

% this function tests the functionality of FT_CONNECTIVITYANALYSIS
% on frequency domain channel data

% apart from using FT_CONNECTIVITYANALYSIS, it also relies on
% FT_CONNECTIVITYSIMULATION, FT_FREQANALYSIS, FT_MVARANALYSIS

global ft_default;
ft_default.feedback = 'no';

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
cfgc.granger.sfmethod  = 'bivariate';
c5b            = ft_connectivityanalysis(cfgc, freq);
cfgc.granger   = rmfield(cfgc.granger, 'sfmethod');
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
cfgc.method    = 'ddtf';
c10            = ft_connectivityanalysis(cfgc, freq);
c10m           = ft_connectivityanalysis(cfgc, mfreq);
cfgc.method    = 'gpdc';
c11            = ft_connectivityanalysis(cfgc, freq);
c11m           = ft_connectivityanalysis(cfgc, mfreq);


cfgc             = [];
cfgc.partchannel = 'signal003'; % this should destroy coherence between 1 and 2
cfgc.method      = 'coh';
c10              = ft_connectivityanalysis(cfgc, freq);

cfgc             = [];
cfgc.method      = 'coh';
cfgc.channelcmb  = {'signal001' 'signal002'};
c11              = ft_connectivityanalysis(cfgc, freq);
cfgc.channelcmb  = {{'signal001'} {'signal002';'signal003'}};
c12              = ft_connectivityanalysis(cfgc, freq);

cfgc             = [];
cfgc.method      = 'granger';
cfgc.channelcmb  = {'signal001' 'signal002'};
c13              = ft_connectivityanalysis(cfgc, freq); %gives a 'chan_chan_freq' matrix
cfgc.sfmethod    = 'bivariate';
c14              = ft_connectivityanalysis(cfgc, freq); %gives a 'chan_freq' matrix
cfgc.channelcmb  = {{'signal001'} {'signal002';'signal003'}};
c15              = ft_connectivityanalysis(cfgc, freq); %gives a 'chan_freq' matrix (4 x nfreq)


% this part tests the functionality of blockwisegranger 
% FIXME as of yet it does not contain any explicit assertions, just the 
% code is there to check whether or not it crashes.
% Checks that can be done are: 
%  - pairwise spectral factorization should yield same results as multivariate when only 2 channels in input
%  - pairwise/multivariate/blockwise should yield same results when only 2 channels and 1 channel per block
%  = multivariate and blockwise should yield same results when 1 channel per block in general

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 channels not connected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first create some data
cfg             = [];
cfg.ntrials     = 500;
cfg.nsignal     = 2;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.method      = 'ar';
cfg.params(:,:,1) = [ 0.8 0; 
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

% connectivityanalysis
cfgc           = [];
cfgc.method    = 'granger';
c1m            = ft_connectivityanalysis(cfgc, mfreq);
c1b            = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.sfmethod = 'bivariate';
c1b2           = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.sfmethod = 'multivariate';
cfgc.granger.block(1).name =  'block1';
cfgc.granger.block(1).label = freq.label(1);
cfgc.granger.block(2).name = 'block2';
cfgc.granger.block(2).label = freq.label(2);
c1xm           = ft_connectivityanalysis(cfgc, mfreq);
c1x            = ft_connectivityanalysis(cfgc, freq);

%--------------------------------------------------------
% now duplicate both channels, and add some noise
data2 = data;
for k = 1:numel(data.trial)
  data2.trial{k} = [1 0;1 0;0 1;0 1]*data.trial{k} + 0.5.*randn(4,size(data.trial{k},2));
end
data2.label = {'signal001a';'signal001b';'signal002a';'signal002b'};

freq2    = ft_freqanalysis(cfgf, data2);
mdata2   = ft_mvaranalysis(cfgm,  data2);
mfreq2   = ft_freqanalysis(cfgfm, mdata2);

cfgc = [];
cfgc.method = 'granger';
c2  = ft_connectivityanalysis(cfgc, freq2);
c2m = ft_connectivityanalysis(cfgc, mfreq2);
cfgc.granger.block(1).name =  'block1';
cfgc.granger.block(1).label = freq2.label(1:2);
cfgc.granger.block(2).name = 'block2';
cfgc.granger.block(2).label = freq2.label(3:4);
c2x = ft_connectivityanalysis(cfgc, freq2);

tmp = c2.grangerspctrm([1 3],[1 3],:) + ...
      c2.grangerspctrm([1 4],[1 4],:) + ...
      c2.grangerspctrm([2 3],[2 3],:) + ...
      c2.grangerspctrm([2 4],[2 4],:);

c2x2 = c2x;
c2x2.grangerspctrm = tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 channels connected 2->1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first create some data
cfg             = [];
cfg.ntrials     = 500;
cfg.nsignal     = 2;
cfg.triallength = 1;
cfg.fsample     = 200;
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
freqsub        = ft_selectdata(freq, 'foilim', freq.freq(2:end));

% connectivityanalysis
cfgc           = [];
cfgc.method    = 'granger';
c1m            = ft_connectivityanalysis(cfgc, mfreq);
c1b            = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.sfmethod = 'bivariate';
c1b2           = ft_connectivityanalysis(cfgc, freq);
c1b2sub        = ft_connectivityanalysis(cfgc, freqsub);
cfgc.granger.sfmethod = 'multivariate';
cfgc.granger.block(1).name =  'block1';
cfgc.granger.block(1).label = freq.label(1);
cfgc.granger.block(2).name = 'block2';
cfgc.granger.block(2).label = freq.label(2);
c1xm           = ft_connectivityanalysis(cfgc, mfreq);
c1x            = ft_connectivityanalysis(cfgc, freq);

%--------------------------------------------------------
% now duplicate both channels, and add some noise
data2 = data;
for k = 1:numel(data.trial)
  data2.trial{k} = [1 0;1 0;0 1;0 1]*data.trial{k} + 0.5.*randn(4,size(data.trial{k},2));
end
data2.label = {'signal001a';'signal001b';'signal002a';'signal002b'};

freq2    = ft_freqanalysis(cfgf, data2);
mdata2   = ft_mvaranalysis(cfgm,  data2);
mfreq2   = ft_freqanalysis(cfgfm, mdata2);


cfgc = [];
cfgc.method = 'granger';
c2  = ft_connectivityanalysis(cfgc, freq2);
c2m = ft_connectivityanalysis(cfgc, mfreq2);
cfgc.granger.block(1).name =  'block1';
cfgc.granger.block(1).label = freq2.label(1:2);
cfgc.granger.block(2).name = 'block2';
cfgc.granger.block(2).label = freq2.label(3:4);
c2x = ft_connectivityanalysis(cfgc, freq2);

tmp = c2.grangerspctrm([1 3],[1 3],:) + ...
      c2.grangerspctrm([1 4],[1 4],:) + ...
      c2.grangerspctrm([2 3],[2 3],:) + ...
      c2.grangerspctrm([2 4],[2 4],:);

c2x2 = c2x;
c2x2.grangerspctrm = tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 channels not connected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create some data
cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';
cfg.params(:,:,1) = [ 0.8 0   0; 
                      0   0.9 0;
                      0   0   0.5];
cfg.params(:,:,2) = [-0.5    0  0; 
                        0 -0.8  0; 
                        0    0 -0.2];
cfg.noisecov      = [0.3 0 0;
                       0 1 0;
                       0 0 0.2];

data            = ft_connectivitysimulation(cfg);

% freqanalysis
cfgf           = [];
cfgf.method    = 'mtmfft';
cfgf.output    = 'fourier';
cfgf.tapsmofrq = 2;
freq           = ft_freqanalysis(cfgf, data);

cfgc             = [];
cfgc.method      = 'granger';
cfgc.granger.sfmethod = 'bivariate';
g1               = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.sfmethod = 'multivariate';
g2               = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.block(1).name = freq.label{1};
cfgc.granger.block(1).label = freq.label(1);
cfgc.granger.block(2).name = freq.label{2};
cfgc.granger.block(2).label = freq.label(2);
cfgc.granger.block(3).name = freq.label{3};
cfgc.granger.block(3).label = freq.label(3);
g3               = ft_connectivityanalysis(cfgc, freq);

cfgc.granger.conditional = 'yes';
g4               = ft_connectivityanalysis(cfgc, freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 channels connected 3->2 and 1->3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create some data
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

% freqanalysis
cfgf           = [];
cfgf.method    = 'mtmfft';
cfgf.output    = 'fourier';
cfgf.tapsmofrq = 2;
freq           = ft_freqanalysis(cfgf, data);

cfgc             = [];
cfgc.method      = 'granger';
cfgc.granger.sfmethod = 'bivariate';
g1               = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.sfmethod = 'multivariate';
g2               = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.block(1).name = freq.label{1};
cfgc.granger.block(1).label = freq.label(1);
cfgc.granger.block(2).name = freq.label{2};
cfgc.granger.block(2).label = freq.label(2);
cfgc.granger.block(3).name = freq.label{3};
cfgc.granger.block(3).label = freq.label(3);
g3               = ft_connectivityanalysis(cfgc, freq);

cfgc.granger.conditional = 'yes';
g4               = ft_connectivityanalysis(cfgc, freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 blocks of channels connected 3->1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create some data
cfg = [];
cfg.method      = 'ar';
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 6;
cfg.bpfilter    = 'no';
cfg.demean       = 'yes';
cfg.params(:,:,1)      = [0.8   0       0     0       0      0   ; 
                          0.1   0.8     0     0       0      0.2   ;
                          0     0       0.9   0       0      0   ;
                          0     0       0     0.9     0      0   ; 
                          0     0       0     0       0.5    0   ;
                          0     0      0     0       0      0.5];
                      
cfg.params(:,:,2)      = [-0.5    0     0     0      0       0   ; 
                           0     -0.5   0     0      0       0   ;
                           0      0    -0.8   0      0       0   ;
                           0      0     0    -0.8    0       0   ;  
                           0      0     0     0     -0.2     0   ;
                           0      0     0     0      0      -0.2];
                       
cfg.noisecov     = [1     0.5     0     0     0     0   ; 
                    0.5   1       0     0     0     0   ;
                    0     0       1     0.5   0     0   ;
                    0     0       0.5   1     0     0   ; 
                    0     0       0     0     1     0.5 ;
                    0     0       0     0     0.5   1  ];
                    
data = ft_connectivitysimulation(cfg);

% freqanalysis
cfgf           = [];
cfgf.method    = 'mtmfft';
cfgf.output    = 'fourier';
cfgf.tapsmofrq = 2;
freq           = ft_freqanalysis(cfgf, data);

cfgc             = [];
cfgc.method      = 'granger';
cfgc.granger.sfmethod = 'bivariate';
g1               = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.sfmethod = 'multivariate';
g2               = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.block(1).name =  'block1';
cfgc.granger.block(1).label = freq.label(1:2);
cfgc.granger.block(2).name = 'block2';
cfgc.granger.block(2).label = freq.label(3:4);
cfgc.granger.block(3).name = 'block3';
cfgc.granger.block(3).label = freq.label(5:6);
g3               = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.conditional = 'yes';
g4               = ft_connectivityanalysis(cfgc, freq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 blocks of channels connected 1->3 and 3->2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create some data
cfg = [];
cfg.method      = 'ar';
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 8;
cfg.bpfilter    = 'no';
cfg.demean       = 'yes';
cfg.params(:,:,1)      = [0.8   0       0     0       0      0      0      0; 
                          0     0.8     0     0       0      0      0      0;
                          0     0       0.9   0       0.5    0.5    0      0;
                          0     0       0     0.9     0.5    0.5    0      0; 
                          0.4   0.4     0     0       0.5    0      0      0;
                          0.4   0.4     0     0       0      0.5    0      0;
                          0     0       0     0       0      0      0.7    0;
                          0     0       0     0       0      0      0      0.7];      
                      
cfg.params(:,:,2)      = [-0.5    0     0     0      0       0      0      0; 
                           0     -0.5   0     0      0       0      0      0;
                           0      0    -0.8   0      0       0      0      0;
                           0      0     0    -0.8    0       0      0      0;  
                           0      0     0     0     -0.2     0      0      0;
                           0      0     0     0      0      -0.2    0      0;
                           0      0     0     0      0       0     -0.4    0;
                           0      0     0     0      0       0      0     -0.4];

cfg.noisecov     = [1     0.5     0     0     0     0    0    0; 
                    0.5   1       0     0     0     0    0    0;
                    0     0       1     0.5   0     0    0    0;
                    0     0       0.5   1     0     0    0    0; 
                    0     0       0     0     1     0.5  0    0;
                    0     0       0     0     0.5   1    0    0;
                    0     0       0     0     0     0    1    0.5;
                    0     0       0     0     0     0    0.5  1];

data = ft_connectivitysimulation(cfg);

% freqanalysis
cfgf           = [];
cfgf.method    = 'mtmfft';
cfgf.output    = 'fourier';
cfgf.tapsmofrq = 2;
freq           = ft_freqanalysis(cfgf, data);

cfgc             = [];
cfgc.method      = 'granger';
cfgc.granger.sfmethod = 'bivariate';
g1               = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.sfmethod = 'multivariate';
g2               = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.block(1).name =  'block1';
cfgc.granger.block(1).label = freq.label(1:2);
cfgc.granger.block(2).name = 'block2';
cfgc.granger.block(2).label = freq.label(3:4);
cfgc.granger.block(3).name = 'block3';
cfgc.granger.block(3).label = freq.label(5:6);
cfgc.granger.block(4).name = 'block4';
cfgc.granger.block(4).label = freq.label(7:8);
g3               = ft_connectivityanalysis(cfgc, freq);
cfgc.granger.conditional = 'yes';
g4               = ft_connectivityanalysis(cfgc, freq);

