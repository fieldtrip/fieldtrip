function test_example_irasa()
% MEM 3gb
% WALLTIME 00:05:00
% DEPENDENCY ft_freqanalysis ft_specest_irasa external/signal/resample.m
% tested on MATLAB_R2022a, macOS_Monterey_12.4, FT_2f387ff (@JM fixed resample scaling issue)

clear all;

% set simulation parameters
A = 1; % scale of 1/f amplitude
C = 1; % 1/f slope
O = 1; % weight of oscillatory components of the simulated data

lf = 1; % lower bound of freq
hf = 300; % higher bound of freq
sl = 600; % spectral lines

fs = 1000; % sampling rate
n = 60000; % time pnts
t = (1:n)/fs; % time axis

% simulate data
for rpt = 1:1

    %     try
    %       % generate pink noise
    %         dspobj = dsp.ColoredNoise('Color', 'pink', 'SamplesPerFrame', length(t));
    %         fn = dspobj()';
    %     catch
    %         % use another method to make pink noise when dsp.ColoredNoise returns licence error
    %         fn = cumsum(randn(1,length(t)));
    %         fn = fn./max(abs(fn)); %%% @JM This scale doesn't seem right?
    %     end

    % another way to simulate pink noise which give users more control over the features of the noise
    freq = linspace(lf, hf, sl); % sampling frequesies
    fn = zeros(size(t));
    for i=1:length(freq) % cummulative sum over freq
        fn = fn + (A * 1/freq(i)^C) * cos(2*pi*freq(i)*t + rand*2*pi); % 1/f amplitude = a*(1/f^c)
    end

    % add a 10Hz and 60 Hz oscillation
    data.trial{1,rpt} = fn + O * cos(2*pi*10*t) + O * cos(2*pi*60*t);
    data.time{1,rpt}  = t;
    data.label{1}     = 'chan';
    data.trialinfo(rpt,1) = rpt;
end

% chunk 2-second segments (gives 1Hz frequency resolution) for long/continous trials
cfg           = [];
cfg.length    = 2; % freqency resolution = 1/2^floor(log2(cfg.length*0.9))
cfg.overlap   = 0.5;
data          = ft_redefinetrial(cfg, data);

% compute the fractal and original spectra
tic
cfg               = [];
cfg.foilim        = [1 200];
cfg.pad           = 'nextpow2';
cfg.method        = 'irasa';
cfg.output        = 'fractal';
fractal = ft_freqanalysis(cfg, data);
cfg.output        = 'original';
original = ft_freqanalysis(cfg, data);
toc % ~28s

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

if A~=0 && O==0
    title('pure fractal signal');
elseif A==0 && O~=0
    title('pure oscillatory signal');
elseif A~=0 && O~=0
    title('mixed signal');
end