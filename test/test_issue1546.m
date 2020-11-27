function test_issue1546(F, O)
% [F O] = weight of (F)ractal and (O)scillatory components of the simulated data

%% simulate data
t = (1:2000)/1000; % time axis
for rpt = 1:100
    % generate pink noise
    dspobj = dsp.ColoredNoise('Color', 'pink', ...
        'SamplesPerFrame', length(t));
    fn = dspobj()';
    
    % add a 10Hz and 60 Hz oscillation
    data.trial{1,rpt} = F * fn + O * cos(2*pi*10*t) + O * cos(2*pi*60*t);
    data.time{1,rpt}  = t;
    data.label{1}     = 'chan';
    data.trialinfo(rpt,1) = rpt;
end

%% test_ft_specset_irasa
cfg               = [];
cfg.foilim        = [1 200];
cfg.pad           = 'nextpow2';
cfg.method        = 'irasa';
cfg.output        = 'fracpow';
frac = ft_freqanalysis(cfg, data);
cfg.output        = 'pow';
orig = ft_freqanalysis(cfg, data);
% display the spectra in log-log scale
figure();
loglog(orig.freq, orig.powspctrm,'r-');hold on;
loglog(frac.freq, frac.powspctrm,'r--');hold on;

% %% Wen&Liu,2016@https://purr.purdue.edu/publications/1987/1
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

%% ft_specset_irasa
% partition the data into ten overlapping sub-segments
w = data.time{1}(end)-data.time{1}(1); % window length
cfg               = [];
cfg.length        = w*.9;
cfg.overlap       = 1-(((w-cfg.length)/cfg.length)/(15-1));%the number of segement = 15;
data = ft_redefinetrial(cfg, data);

% perform IRASA and regular spectral analysis
cfg               = [];
cfg.foilim        = [1 200];
cfg.taper         = 'hanning';
cfg.pad           = 'nextpow2';
cfg.method        = 'irasa';
frac = ft_freqanalysis_old(cfg, data);
cfg.method        = 'mtmfft';
orig = ft_freqanalysis_old(cfg, data);
% display the spectra in log-log scale
loglog(orig.freq, orig.powspctrm,'b-');hold on;
loglog(frac.freq, frac.powspctrm,'b--');hold on;

xlabel('log-freq'); ylabel('log-power');
% legend({'orgi test\_ft\_irasa','frac test\_ft\_irasa','orgi ft\_irasa','frac ft\_irasa','orgi wen','frac wen'},'location','southeast');
legend({'orgi test\_ft\_irasa','frac test\_ft\_irasa','orgi ft\_irasa','frac ft\_mtmfft'},'location','southwest');
if F~=0&&O==0
    title('pure fractal signal');
elseif F==0&&O~=0
    title('pure oscillatory signal');
elseif F~=0&&O~=0
    title('mixed signal');
end

