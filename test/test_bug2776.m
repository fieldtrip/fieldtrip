function test_bug2776

% WALLTIME 00:10:00
% MEM 1gb

% see http://nl.mathworks.com/help/signal/ug/psd-estimate-using-fft.html

Fs = 1000;
t = (1:Fs)/Fs;
% x = randn(size(time));
x = cos(2*pi*100*t) + randn(size(t));

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

%% fieldtrip style

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

