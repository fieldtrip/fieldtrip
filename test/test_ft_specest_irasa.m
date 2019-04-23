function test_ft_specest_irasa

% MEM 1500mb
% WALLTIME 00:10:00

% DEPENDENCY ft_freqanalysis

% script demonstrating the extraction of rhythmic spectral 
% features from the electrophysiological signal


% generate trials with a 15 Hz oscillation embedded in pink noise
t = (1:1:1000)/1000; % time axis
for rpt = 1:1000
  % generate pink noise
  dspobj = dsp.ColoredNoise('Color', 'pink', ...
    'SamplesPerFrame', length(t));
  fn = dspobj()';
  
  % add a 15 Hz oscillation
  data.trial{1,rpt} = fn + cos(2*pi*15*t); 
  data.time{1,rpt} = t;
  data.label{1} = 'chan';
end

% perform regular frequency analysis and IRASA
cfg = [];
cfg.foilim = [1 50];
cfg.taper = 'hanning';
cfg.method = 'irasa';
frac = ft_freqanalysis(cfg, data);
cfg.method = 'mtmfft';
orig = ft_freqanalysis(cfg, data);

% subtract the fractal component from the power spectrum
cfg = [];
cfg.parameter = 'powspctrm';
cfg.operation = 'x2-x1';
osci = ft_math(cfg, frac, orig);

% plot the fractal component and the power spectrum 
figure; plot(frac.freq, frac.powspctrm, ...
  'linewidth', 3, 'color', [0 0 0])
hold on; plot(orig.freq, orig.powspctrm, ...
  'linewidth', 3, 'color', [.6 .6 .6])

% plot the full-width half-maximum of the oscillatory component
f = fit(osci.freq', osci.powspctrm', 'gauss1');
mean = f.b1;
std = f.c1/sqrt(2)*2.3548;
fwhm = [mean-std/2 mean+std/2];
yl = get(gca, 'YLim');
p = patch([fwhm flip(fwhm)], [yl(1) yl(1) yl(2) yl(2)], [1 1 1]);
uistack(p, 'bottom');
legend('FWHM oscillation', 'Fractal component', 'Power spectrum');
xlabel('Frequency'); ylabel('Power');
