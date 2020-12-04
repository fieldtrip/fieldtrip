function test_issue1546

% MEM 3gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqanalysis ft_specest_irasa ft_specest_irasa_new

%% simulate data
F = 1; % weight of (F)ractal components of the simulated data
O = 1; % weight of (O)scillatory components of the simulated data

t = (1:60000)/1000; % time axis
for rpt = 1:1
    try
      % generate pink noise
        dspobj = dsp.ColoredNoise('Color', 'pink', 'SamplesPerFrame', length(t));
        fn = dspobj()';
    catch
        % use another method to make pink noise when dsp.ColoredNoise returns licence error
        fn = cumsum(randn(1,length(t))); 
    end
    
    % add a 10Hz and 60 Hz oscillation
    data.trial{1,rpt} = F * fn + O * cos(2*pi*10*t) + O * cos(2*pi*60*t);
    data.time{1,rpt}  = t;
    data.label{1}     = 'chan';
    data.trialinfo(rpt,1) = rpt;
end

%% chunk 2-second segments (gives 1Hz frequency resolution) for long/continous trials
cfg           = [];
cfg.length    = 2; % freqency resolution = 1/2^floor(log2(cfg.length*0.9))
cfg.overlap   = 0.5;
data          = ft_redefinetrial(cfg, data);

%% ft_specest_irasa_new
cfg               = [];
cfg.foilim        = [1 200];
cfg.pad           = 'nextpow2';
cfg.method        = 'irasa_new';
cfg.output        = 'fractal';
fractal = ft_freqanalysis(cfg, data);
cfg.output        = 'mixed';
original = ft_freqanalysis(cfg, data);

% display the spectra in log-log scale
figure();
 hold on;
loglog(original.freq, original.powspctrm, 'r-');
loglog(fractal.freq, fractal.powspctrm, 'r--');

%% Wen&Liu,2016@https://purr.purdue.edu/publications/1987/1
% % convert data structure
% dat = [];
% for i=data.trialinfo'
%     dat = [dat; data.trial{i}];
% end
% dat = dat';
% 
% % set parameter
% srate = 1000; % sampling frequency
% frange = [1 200];
% 
% % separate fractal and oscillatory components
% frac = amri_sig_fractal(dat,srate,'frange',frange);
% 
% % display the spectra in log-log scale
% loglog(frac.freq, mean(frac.mixd,2),'g-');hold on;
% loglog(frac.freq, mean(frac.frac,2),'g--');hold on;

%% Use the implementation that was made by Arjen
% partition the data into ten overlapping sub-segments
w = data.time{1}(end)-data.time{1}(1); % window length
cfg               = [];
cfg.length        = w*.9;
cfg.overlap       = 1-(((w-cfg.length)/cfg.length)/(10-1)); % the number of segement = 10;
data = ft_redefinetrial(cfg, data);

% perform IRASA and regular spectral analysis
cfg               = [];
cfg.foilim        = [1 200];
cfg.taper         = 'hanning';
cfg.pad           = 'nextpow2';
cfg.method        = 'irasa';
frac = ft_freqanalysis(cfg, data);
cfg.method        = 'mtmfft';
original = ft_freqanalysis(cfg, data);

% display the spectra in log-log scale
loglog(original.freq, original.powspctrm, 'b-');
loglog(fractal.freq, fractal.powspctrm, 'b--');

xlabel('log-freq'); ylabel('log-power');
legend({'mixed new\_ft\_irasa','frac new\_ft\_irasa','mixed ft\_mtmfft','frac ft\_irasa'},'location','southwest');

if F~=0 && O==0
    title('pure fractal signal');
elseif F==0 && O~=0
    title('pure oscillatory signal');
elseif F~=0 && O~=0
    title('mixed signal');
end
