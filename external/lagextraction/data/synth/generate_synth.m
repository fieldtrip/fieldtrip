function [EEG,ref,latency,snr] = generate_synth(EEG,options)
%   GENERATE_SYNTH   
%       [EEG,REF,LATENCY,SNR] = GENERATE_SYNTH(EEG,OPTIONS)
% 
%   Created by Alexandre Gramfort on 2009-05-10.
%   Copyright (c)  Alexandre Gramfort. All rights reserved.

% $Id: generate_synth.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

if nargin<2
    options.null = 0;
end

if ~isfield(options, 'mean_latency')
    options.mean_latency = 0.3;
end
mean_latency = options.mean_latency;

if ~isfield(options, 'mean_sigma')
    options.mean_sigma = 0.05;
end
mean_sigma = options.mean_sigma;

if ~isfield(options, 'mean_frequency')
    options.mean_frequency = 5;
end
mean_frequency = options.mean_frequency;

if ~isfield(options, 'mean_strech')
    options.mean_strech = 0;
end
mean_strech = options.mean_strech;

if ~isfield(options, 'std_latency')
    options.std_latency = 0;
end
std_latency = options.std_latency;

if ~isfield(options, 'std_sigma')
    options.std_sigma = 0;
end
std_sigma = options.std_sigma;

if ~isfield(options, 'std_frequency')
    options.std_frequency = 0;
end
std_frequency = options.std_frequency;

if ~isfield(options, 'std_strech')
    options.std_strech = 0;
end
std_strech = options.std_strech;

if ~isfield(options, 'snr')
    options.snr = 2;
end
snr = options.snr;

if ~isfield(options, 'use_ar')
    options.use_ar = false;
end
use_ar = options.use_ar;

if isfield(options, 'null')
    options = rmfield(options,'null');
end

nepochs = EEG.trials;

R = TruncatedGaussian(1,[-2 2],[nepochs,1]);

latency   = mean_latency   + std_latency * R;
sigma     = mean_sigma     + std_sigma * R;
frequency = mean_frequency + std_frequency * R;
strech    = mean_strech    + std_strech * R;

% fname_suffix = [ '_', ...
%                  'mfreq_',num2str(mean_frequency),'_', ...
%                  'stdfreq_',num2str(std_frequency),'_', ...
%                  'mlat_',num2str(mean_latency),'_', ...
%                  'stdlat_',num2str(std_latency),'_', ...
%                  'msigma_',num2str(mean_sigma),'_', ...
%                  'stdsigma_',num2str(std_sigma),'_', ...
% ];

% noise AR filter config
b = fir1(1024, .5);
[d,p0] = lpc(b,7);
noise = zeros([1 EEG.pnts EEG.trials]);

dphase = pi/3;
gabor = @(time,freq,sigma,dphase) cos(2*pi*freq*time - dphase) .* exp(-time.^2./2./sigma.^2);
ref = gabor(EEG.times./ (1 + mean_strech)-mean_latency,mean_frequency,mean_sigma,dphase);
ref = ref ./ norm(ref);

for iepoch = 1:nepochs
    time_epoch = EEG.times ./ (1 + strech(iepoch)) - latency(iepoch);
    EEG.data(1,:,iepoch) = gabor(time_epoch,frequency(iepoch),sigma(iepoch),dphase);
    EEG.data(1,:,iepoch) = EEG.data(1,:,iepoch) ./ norm(EEG.data(1,:,iepoch));

    if use_ar
        noise(1,:,iepoch) =  filter(1,d,sqrt(p0)*randn(EEG.pnts,1))'; % AR noise
    else
        noise(1,:,iepoch) =  randn(1,size(noise,2)); % white noise
    end
    noise(1,:,iepoch) = noise(1,:,iepoch) ./ norm(noise(1,:,iepoch));
end

noise = noise ./ norm(noise(:));
noise = norm(EEG.data(:)) * noise ./ snr;

EEG.data = EEG.data + noise;

EEG.times = EEG.times;
EEG.setname = 'synth_realignment';
EEG.xmin = 0;
EEG.xmax = 1;
EEG.icaact = [];
EEG.icawinv = [];
EEG.icasphere = [];
EEG.icaweights = [];

end %  function