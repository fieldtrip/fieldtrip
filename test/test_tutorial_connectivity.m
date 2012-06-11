function test_tutorial_connectivity

% TEST test_tutorial_connectivity
% TEST ft_connectivityanalysis ft_connectivitysimulation ft_freqanalysis ft_connectivityplot ft_mvaranalysis

% disable verbose output
global ft_default;
ft_default.feedback = 'no';

%% simulate data
cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';
cfg.params(:,:,1) = [ 0.8    0    0 ;
  0  0.9  0.5 ;
  0.4    0  0.5];

cfg.params(:,:,2) = [-0.5    0    0 ;
  0 -0.8    0 ;
  0    0 -0.2];

cfg.noisecov      = [ 0.3    0    0 ;
  0    1    0 ;
  0    0  0.2];
data              = ft_connectivitysimulation(cfg);


%% mvaranalysis
cfg         = [];
cfg.order   = 5;
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, data);

%% freqanalysis 1
cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);

%% freqanalysis 2
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.tapsmofrq = 2;
freq          = ft_freqanalysis(cfg, data);

%% connectivityanalysis
cfg           = [];
cfg.method    = 'coh';
coh           = ft_connectivityanalysis(cfg, freq);
cohm          = ft_connectivityanalysis(cfg, mfreq);

%% visualisation
cfg           = [];
cfg.parameter = 'cohspctrm';
ft_connectivityplot(cfg, coh, cohm);

%% do the same for granger
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, mfreq);

cfg           = [];
cfg.parameter = 'grangerspctrm';
ft_connectivityplot(cfg, granger);

figure
for row=1:3
  for col=1:3
    subplot(3,3,(row-1)*3+col);
    plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)))
    ylim([0 1])
  end
end


