function test_ft_specest_irasa

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_freqanalysis ft_specest_irasa

% We do not have that many concurrent licenses of the DSP Systems toolbox (https://nl.mathworks.com/products/dsp-system.html)
% which was used in the initial version of this script. Therefore this test script does not use the DSP toolbox by default.
usedsptoolbox = false;

% simulate data
t = (1:2000)/1000; % time axis
for rpt = 1:20
    if usedsptoolbox
      % generate pink noise using the DSP toolbox
        dspobj = dsp.ColoredNoise('Color', 'pink', 'SamplesPerFrame', length(t));
        fn = dspobj()';
    else
        % use another method to make pink noise
        fn = cumsum(randn(1,length(t))); 
    end
    
    % add line noise
    data.trial{1,rpt} = fn + cos(2*pi*50*t) + cos(2*pi*100*t);
    data.time{1,rpt}  = t;
    data.label{1}     = 'chan';
    data.trialinfo(rpt,1) = rpt;
end

% using unfiltered data
cfg = [];
cfg.method = 'irasa';
cfg.output = 'original';
cfg.pad   = 'nextpow2';
freq = ft_freqanalysis(cfg, data);
cfg.output = 'fractal';
freqI = ft_freqanalysis(cfg, data);

% what happens if we use bandpassfiltered data?
cfg2            = [];
cfg2.bpfilter   = 'yes';
cfg2.bpfilttype = 'firws';
cfg2.bpfreq     = [60 150];
datafilt       = ft_preprocessing(cfg2, data);

cfg2.bpfilter = 'no';
cfg2.hpfilter = 'yes';
cfg2.hpfreq   = 60;
cfg2.hpfilttype = 'firws';
datafilt2       = ft_preprocessing(cfg2, data);

cfg2.bpfilter   = 'no';
cfg2.hpfilter   = 'no';
cfg2.dftfilter  = 'yes'; % notch filter to filter out line noise
cfg2.dftfreq    = [50 100];
datafilt3       = ft_preprocessing(cfg2, data);

cfg.output = 'original';
freqfilt = ft_freqanalysis(cfg, datafilt);
freqfilt2 = ft_freqanalysis(cfg, datafilt2);
freqfilt3 = ft_freqanalysis(cfg, datafilt3);

cfg.output = 'fractal';
freqfiltI = ft_freqanalysis(cfg, datafilt);
freqfiltI2 = ft_freqanalysis(cfg, datafilt2);
freqfiltI3 = ft_freqanalysis(cfg, datafilt3);

figure; 
semilogy(freq.freq, [freq.powspctrm;freqfilt.powspctrm;freqfilt2.powspctrm;freqfilt3.powspctrm;...
                     freqI.powspctrm;freqfiltI.powspctrm;freqfiltI2.powspctrm;freqfiltI3.powspctrm]);
legend({'orig-unfilterd';'orig-bpfiltered';'orig-hpfiltered';'orig-dftfiltered';...
        'frac-unfiltered';'frac-bpfiltered';'frac-hpfiltered';'frac-dftfiltered'});

