function [EEG,ref,latencies] = generate_synth_eeg_data(EEG,options)

% $Id: generate_synth_eeg_data.m 2 2009-06-16 19:24:10Z gramfort $
% $LastChangedBy: gramfort $
% $LastChangedDate: 2009-06-16 15:24:10 -0400 (Mar, 16 jui 2009) $
% $Revision: 2 $

if nargin<2
    options.null = 0;
end

if ~isfield(options, 'snr')
    options.snr = 2;
end
snr = options.snr;

if ~isfield(options, 'sigma_latencies')
    options.sigma_latencies = 0.05;
end
sigma_latencies = options.sigma_latencies;

if ~isfield(options, 't1')
    options.t1 = 0.15;
end
t1 = options.t1;

if ~isfield(options, 'use_randn')
    options.use_randn = true;
end
use_randn = options.use_randn;

if ~isfield(options, 'use_ar')
    options.use_ar = false;
end
use_ar = options.use_ar;

if isfield(options, 'null')
    options = rmfield(options,'null');
end

tmin = 1;
tmax = EEG.pnts;

% Shape of the response
d1 = 0.05;
m1 = 1;
ref = (EEG.times-t1(1)) .* exp(-(EEG.times-t1(1))/d1) .* (sign(EEG.times-t1(1))+1)/2 / d1 * exp(1) * m1;
ref = ref ./ norm(ref);

EEG.data = zeros([1 EEG.pnts EEG.trials]);
noise = zeros([1 EEG.pnts EEG.trials]);

% noise AR filter config
b = fir1(1024, .5);
[d,p0] = lpc(b,7);

if use_randn
    latencies = sigma_latencies*randn([1,EEG.trials]);
else
    latencies = 2*sigma_latencies*(rand([1,EEG.trials])-0.5);
end

for k=1:EEG.trials
    deltat = t1(round(unifrnd(1,length(t1))));
    % deltat = t1(randint(1,1,[1,length(t1)]));
    latencies(k) = latencies(k) + deltat;
    lag = latencies(k);
    EEG.data(1,:,k) = (EEG.times-lag) .* exp(-(EEG.times-lag)/d1) .* ...
                      (sign(EEG.times-lag)+1)/2 / d1 * exp(1) * m1;
    EEG.data(1,:,k) = EEG.data(1,:,k) ./ norm(EEG.data(1,:,k));
    if use_ar
        noise(1,:,k) =  filter(1,d,sqrt(p0)*randn(EEG.pnts,1))'; % AR noise
    else
        noise(1,:,k) =  randn(1,size(noise,2)); % white noise
    end
    noise(1,:,k) = noise(1,:,k) ./ norm(noise(1,:,k));
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
