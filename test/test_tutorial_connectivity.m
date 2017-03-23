function test_tutorial_connectivity

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_connectivityanalysis ft_connectivitysimulation ft_freqanalysis ft_connectivityplot ft_mvaranalysis

% This is the first section of the connectivity tutorial, which
% starts with an MVAR model and then uses parametric and nonparametric 
% spectral decomposition for coherence and granger

% See also test_tutorial_connectivity2 and test_tutorial_connectivity3

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

% % add instantenous 50Hz noise
% for tr=1:numel(data.trial)
%   r = rand;
% %   noise = sin((1+r)*pi*50.*data.time{tr}) + cos((1-r)*pi*50.*data.time{tr})
%   for c=1:cfg.nsignal
%     timeshift = -0.002; %in seconds
% %     phaseshift = c/pi/5;
%     phaseshift = 0;
% %     phaseshift = 2*pi/(1/50/(c*timeshift));
%     noise(c, :) = sin(2*pi*(50).*data.time{tr}+phaseshift);
%     for h=[-5 -4 -3 -2 -1 1 2 3 4 5]
% %       phaseshift = 2*pi/(1/(50-h)/(c*timeshift));
% %       noise(c, :) = noise(c, :)+sin(1-r*pi.*(50-h).*data.time{tr}+phaseshift);
%     end
%   end
% %   noise = noise*max(data.trial{tr}(:));
% %   data.trial{tr} = data.trial{tr} + repmat(noise, cfg.nsignal, 1);
%   data.trial{tr} = data.trial{tr} + noise;
% end

% cfg = [];
% cfg.dftfilter =  'yes';
% cfg.dftfreq   = [49.4:0.2:50.6 99.4:.2:100.6 149.4:.2:150.6];
% data = ft_preprocessing(cfg, data);
% cfg.bsfilter =  'yes';
% cfg.bsfreq   = [45 50];
% cfg.bpfilttype = 'fir';
% data = ft_preprocessing(cfg, data);
%   
figure
plot(data.time{1}, data.trial{1}) 
legend(data.label)
xlabel('time (s)')

cfg = [];
cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
ft_databrowser(cfg, data);



% %% resample
% cfg = [];
% cfg.lpfilter = 'yes';
% cfg.lpfreq   = 30;
% data = ft_preprocessing(cfg, data);
% 
% cfg = [];
% cfg.resamplefs = 100;
% cfg.detrend = 'yes';
% data = ft_resampledata(cfg, data);

%% mvaranalysis
cfg         = [];
cfg.order   = 10;
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
cfg.pad       = 4;
freq          = ft_freqanalysis(cfg, data);


figure;
cfg = [];
cfg.channel = 'signal001';
ft_singleplotER(cfg, ft_checkdata(freq, 'cmbrepresentation', 'sparsewithpow'))
%% connectivityanalysis
cfg           = [];
cfg.method    = 'coh';
coh           = ft_connectivityanalysis(cfg, freq);
cohm          = ft_connectivityanalysis(cfg, mfreq);

%% visualisation
cfg           = [];
cfg.parameter = 'cohspctrm';
figure;
ft_connectivityplot(cfg, coh, cohm);

%% do the same for granger
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq);
mgranger       = ft_connectivityanalysis(cfg, mfreq);

cfg           = [];
cfg.parameter = 'grangerspctrm';
figure
ft_connectivityplot(cfg, granger);
figure
for row=1:3
  for col=1:3
    subplot(3,3,(row-1)*3+col);
    plot(granger.freq, smooth(squeeze(granger.grangerspctrm(row,col,:))))
    hold on;
    plot(mgranger.freq, squeeze(mgranger.grangerspctrm(row,col,:)), 'r')
    hold off;
    ylim([0 1])
  end
end


