function test_example_irasa()

% MEM 4gb
% WALLTIME 00:10:00
% DEPENDENCY ft_specest_irasa ft_freqanalysis

% simulate data
F = 1; % weight of (F)ractal components of the simulated data
O = 1; % weight of (O)scillatory components of the simulated data
t = (1:60000)/1000; % time axis
for rpt = 1:1
  % use a simple method to make pink noise
  % that does not rely on the digital signal processing toolbox
  fn = cumsum(randn(1,length(t)));
  fn = fn./max(abs(fn));
  
  % add a 10 Hz and 60 Hz oscillation
  data.trial{1,rpt} = F * fn + O * cos(2*pi*10*t) + O * cos(2*pi*60*t);
  data.time{1,rpt}  = t;
  data.label{1}     = 'chan';
  data.trialinfo(rpt,1) = rpt;
end

% chunk into 2-second segments
cfg               = [];
cfg.length        = 2;
cfg.overlap       = 0.5;
data              = ft_redefinetrial(cfg, data);

% compute the fractal and original spectra
cfg               = [];
cfg.foilim        = [1 200];
cfg.pad           = 'nextpow2';
cfg.method        = 'irasa';
cfg.output        = 'fractal';
fractal = ft_freqanalysis(cfg, data);
cfg.output        = 'original';
original = ft_freqanalysis(cfg, data);

% subtract the fractal component from the power spectrum
cfg               = [];
cfg.parameter     = 'powspctrm';
cfg.operation     = 'x2-x1';
oscillatory = ft_math(cfg, fractal, original);

% display the spectra in log-log scale
figure();
hold on;
plot(log(original.freq), log(original.powspctrm),'k');
plot(log(fractal.freq), log(fractal.powspctrm));
plot(log(fractal.freq), log(oscillatory.powspctrm));
xlabel('log-freq'); ylabel('log-power');
legend({'original','fractal','oscillatory'},'location','southwest');
if F~=0 && O==0
  title('pure fractal signal');
elseif F==0 && O~=0
  title('pure oscillatory signal');
elseif F~=0 && O~=0
  title('mixed signal');
end


% 
% % load some real raw data
% load('S5_raw_segmented.mat')
% 
% % filter and re-reference the raw data
% cfg               = [];
% cfg.hpfilter      = 'yes';
% cfg.hpfiltord     = 3;
% cfg.hpfreq        = 1;
% cfg.lpfilter      = 'yes';
% cfg.lpfiltord     = 3;
% cfg.lpfreq        = 249;
% cfg.bsfilter      = 'yes';
% cfg.bsfiltord     = 3;
% cfg.bsfreq        = [49 51; 99 101; 149 151; 199 201]; % EU line noise
% cfg.reref         = 'yes';
% cfg.refchannel    = 'all';
% data_filt = ft_preprocessing(cfg, data);
% 
% % select electrodes placed over the sensorimotor cortex
% cfg               = [];
% cfg.channel       = {'chan007','chan008','chan015','chan016','chan024', ...
%   'chan092','chan093','chan094','chan095', ...
%   'chan113','chan115','chan117','chan118','chan119'};
% data_sel = ft_selectdata(cfg, data_filt);
% %
% %* Run IRASA with sliding window for the continuous recording
% %
% % segment the data into one-second overlapping chunks
% cfg               = [];
% cfg.length        = 1;
% cfg.overlap       = .95; % 95 percent overlap (sliding 50 ms)
% data_redef = ft_redefinetrial(cfg, data_sel);
% 
% % compute the fractal and original spectra
% cfg               = [];
% cfg.foilim        = [1 50];
% cfg.pad           = 'nextpow2';
% cfg.method        = 'irasa';
% cfg.output        = 'fractal';
% fractal = ft_freqanalysis(cfg, data_redef);
% cfg.output        = 'original';
% original = ft_freqanalysis(cfg, data_redef);
% 
% % subtract the fractal component from the power spectrum
% cfg               = [];
% cfg.parameter     = 'powspctrm';
% cfg.operation     = 'x2-x1';
% oscillatory = ft_math(cfg, fractal, original);
% 
% % extract alpha and beta frequency bands
% figure;
% plot(oscillatory.freq, mean(oscillatory.powspctrm), ...
%   'linewidth', 3, 'color', [.3 .3 .3])
% f = fit(oscillatory.freq', mean(oscillatory.powspctrm)', 'gauss3');
% alpha_band = [f.b1-2 f.b1+2];
% beta_band  = [f.b3-3 f.b3+3];
% yl = get(gca, 'YLim');
% p1 = patch([alpha_band flip(alpha_band)], [yl(1) yl(1) yl(2) yl(2)], [.9 .9 .9]);
% p2 = patch([beta_band flip(beta_band)], [yl(1) yl(1) yl(2) yl(2)], [.8 .8 .8]);
% uistack(p2, 'bottom'); uistack(p1, 'bottom');
% legend('alpha band', 'beta band', 'oscillatory component');
% xlabel('Frequency'); ylabel('Power');
% set(gca, 'YLim', yl);
% %
% %
% %
% %* Localizing spectral features in the sensorimotor cortex
% %
% % read in the cortical surface
% cortex = ft_read_headshape('S5_lh.pial');
% 
% % plot the spatial distribution of alpha rhythmic activity
% cfg               = [];
% cfg.frequency     = alpha_band;
% cfg.avgoverfreq   = 'yes';
% alpha = ft_selectdata(cfg, oscillatory);
% 
% cfg               = [];
% cfg.funparameter  = 'powspctrm';
% cfg.funcolorlim   = 'zeromax';
% cfg.method        = 'surface';
% cfg.interpmethod  = 'sphere_weighteddistance';
% cfg.sphereradius  = 10;
% cfg.camlight      = 'no';
% cfg.funcolormap   = 'parula';
% cfg.colorbar      = 'no';
% ft_sourceplot(cfg, alpha, cortex);
% view([-90 20]);
% material dull; lighting gouraud; camlight
% 
% ft_plot_sens(alpha.elec, 'elecshape', 'disc', 'facecolor', [0 0 0])
% %
% %
% % plot the spatial distribution of beta rhythmic activity
% cfg               = [];
% cfg.frequency     = beta_band;
% cfg.avgoverfreq   = 'yes';
% beta = ft_selectdata(cfg, oscillatory);
% 
% cfg               = [];
% cfg.funparameter  = 'powspctrm';
% cfg.funcolorlim   = 'zeromax';
% cfg.method        = 'surface';
% cfg.interpmethod  = 'sphere_weighteddistance';
% cfg.sphereradius  = 10;
% cfg.camlight      = 'no';
% cfg.funcolormap   = 'parula';
% cfg.colorbar      = 'no';
% ft_sourceplot(cfg, beta, cortex);
% view([-90 20]);
% material dull; lighting gouraud; camlight
% 
% ft_plot_sens(beta.elec, 'elecshape', 'disc', 'facecolor', [0 0 0])
