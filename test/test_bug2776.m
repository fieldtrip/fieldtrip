function test_bug2776

% WALLTIME 00:10:00
% MEM 1gb

% see http://nl.mathworks.com/help/signal/ug/psd-estimate-using-fft.html

Fs = 1000;
t1 = (1:Fs)/Fs;
x1 = cos(2*pi*100*t1) + randn(size(t1));

Fs = 1000;
t2 = (1:(2*Fs))/Fs;
x2 = cos(2*pi*100*t2) + randn(size(t2));

t = t1;
x = x1;

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;

figure
plot(t, x);
title('Original signal')
xlabel('Time (s)')
ylabel('amplitude (au)')

figure
plot(freq,10*log10(psdx))
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')

figure
plot(freq,psdx)
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (au^2/Hz)')

figure
periodogram(x,rectwin(length(x)),length(x))

%% FieldTrip style

data = [];
data.time{1}  = t;
data.trial{1} = x;
data.label    = {'A1'};

cfg = [];
cfg.taper = 'boxcar';
cfg.method = 'mtmfft';
cfg.polyremoval = -1; % no removal of baseline or drift
freq = ft_freqanalysis(cfg, data);

figure
plot(freq.freq,freq.powspctrm)
grid on
title('fieldtrip freqanalysis')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (au^2/Hz)')

assert(all((freq.powspctrm-psdx)<10*eps))

%% let us now explore the effect of padding

cfg = [];
cfg.taper = 'boxcar';
cfg.method = 'mtmfft';
cfg.polyremoval = -1; % no removal of baseline or drift
cfg.pad = 1;
freq1 = ft_freqanalysis(cfg, data);
cfg.pad = 2;
freq2 = ft_freqanalysis(cfg, data);
cfg.pad = 3;
freq3 = ft_freqanalysis(cfg, data);
figure
plot(freq1.freq, freq1.powspctrm, 'b.-'); hold on
plot(freq2.freq, freq2.powspctrm, 'r.-'); hold on
plot(freq3.freq, freq3.powspctrm, 'g.-'); hold on

% the total power should remain the same
assert(abs(sum(freq1.powspctrm)-sum(freq2.powspctrm))<100*eps);
assert(abs(sum(freq2.powspctrm)-sum(freq3.powspctrm))<100*eps);

%% compare narrow and broadband

x1 = cos(2*pi*100*t1) + randn(size(t1));
x2 = [x1 x1];

data_sharp = [];
data_sharp.time{1}  = t1; % 1 second
data_sharp.time{2}  = t2; % 2 seconds
data_sharp.trial{1} = x1;
data_sharp.trial{2} = x2;
data_sharp.label    = {'A1'};

x1 = ft_preproc_bandpassfilter(randn(size(t1)), Fs, [90 110]);
x2 = [x1 x1];

data_broad = [];
data_broad.time{1}  = t1; % 1 second
data_broad.time{2}  = t2; % 2 seconds
data_broad.trial{1} = x1;
data_broad.trial{2} = x2;
data_broad.label    = {'A1'};

cfg = [];
cfg.taper = 'boxcar';
% cfg.taper = 'hanning';
% cfg.taper = 'dpss';
% cfg.tapsmofrq = 5;
cfg.method = 'mtmfft';
cfg.polyremoval = -1; % no removal of baseline or drift
cfg.keeptrials = 'yes';
freq_sharp = ft_freqanalysis(cfg, data_sharp);
freq_broad = ft_freqanalysis(cfg, data_broad);

figure
plot(freq_sharp.freq, squeeze(freq_sharp.powspctrm), '.-')
xlim([80 120])
figure
plot(freq_broad.freq, squeeze(freq_broad.powspctrm), '.-')
xlim([80 120])

assert(diff(sum(freq_sharp.powspctrm,3))<100*eps)
assert(diff(sum(freq_broad.powspctrm,3))<100*eps)
