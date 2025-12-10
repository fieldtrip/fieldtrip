function NSxPowerSpectrum(NSx, channelNumber, colorCode)

% NSxPowerSpectrum
%
% Plots a power spectrum of a given channel.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Use NSxPowerSpectrum
%
% Use settingsManager(NSx, channelNumber, colorCode)
%
%   NSx:           The NSx file containing the data.
%
%   channelNumber: The channel to be plotted.
%
%   colorCode:     The color of the plot. This follows the standard "plot"
%                  colors in MATLAB. See "help plot" for more information.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Kian Torab
%   support@blackrockmicro.com
%   Blackrock Microsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version History
%
% 1.1.0.0: October 2, 2020
%   - Fixed a bug where the sampling frequency is now read from the header
%     file instead of it being fixed at 300 Hz.
%
% 1.1.1.0: October 23, 2020
%   - Fixed a small bug with double defining the function name.
%


if ~exist('colorCode', 'var')
    colorCode = 'b';
end

x = double(NSx.Data(channelNumber,:));

rng default;
Fs = NSx.MetaTags.SamplingFreq;
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