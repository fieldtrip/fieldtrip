function NSxPowerSpectrum(NSx, channelNumber, colorCode)

if ~exist('colorCode', 'var')
    colorCode = 'b';
end

x = double(NSx.Data(channelNumber,:));

rng default;
Fs = 300;
t = linspace(0,1,length(x));

N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)).*abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;
plot(freq,10*log10(psdx), colorCode); grid on;
title('Periodogram Using FFT');
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (dB/Hz)');