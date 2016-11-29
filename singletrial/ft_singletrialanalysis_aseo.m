function [output] = ft_singletrialanalysis_aseo(cfg, data, erp_fft)

% FT_SINGLETRIALANALYSIS_ASEO executes single-trial analysis, using the 
% ASEO algorithm (Xu et al, 2009)
%
% Use as:
% [output] = ft_singletrialanalysis_aseo(cfg, data_fft, erp_fft)
% where data_fft is the observed data in the frequency domain, erp_fft
% contains the initial ERP components in the frequency domain. cfg is a
% configuration structure according to
%
% OUTPUT----
% amp_est    : Estimates of ERP components' amplitude
% lat_est    : Estimates of ERP components' latency
% erp_est    : Estimates of ERP waveforms in time domain
% ar         : Estimated AR coefficients of on-going activity
% noise      : Power spectrum of on-going activity fitted in AR model
% sigma      : Power of the input white noise of AR model for on-going activity
% residual   : Residual signal after removing ERPs in time domain
% rejectflag : Each element of rejectflag indicating that the corresponding
%              trial should be rejected or not. For example, rejectflag(9)==1 means 
%              the 9th trial is rejected.
% corr_est    : Correlation between the original data and the recovered signal

% define some options locally
maxOrderAR      = ft_getopt(cfg.aseo, 'maxOrderAR', 10);
fsample         = ft_getopt(cfg.aseo, 'fsample');
period          = ft_getopt(cfg.aseo, 'sampPeri');
jitter          = ft_getopt(cfg.aseo, 'jitter'); % Latency search window defined in main_ASEO.m
searchGrid      = ft_getopt(cfg.aseo, 'searchGrid'); % Latency search step
thresholdAmpH   = ft_getopt(cfg.aseo, 'thresholdAmpH'); % maximum acceptable amplitude of trials, times of avergae amplitude
thresholdAmpL   = ft_getopt(cfg.aseo, 'thresholdAmpL'); % minimum  acceptable amplitude of trials, times of avergae amplitude
thresholdCorr   = ft_getopt(cfg.aseo, 'thresholdCorr'); % minimum correlation with the original data
nsample         = ft_getopt(cfg.aseo, 'nsmp'); % the number of original samples in the time domain
amp_est         = ft_getopt(cfg.aseo, 'amp_est', []);
erp_est         = ft_getopt(cfg.aseo, 'erp_est', []);
lat_est         = ft_getopt(cfg.aseo, 'lat_est', []);
noise           = ft_getopt(cfg.aseo, 'noise', []);
numiteration    = ft_getopt(cfg.aseo, 'numiteration', 1);
ntrl            = ft_getopt(cfg.aseo, 'ntrl', size(data, 2));
tapsmofrq       = ft_getopt(cfg.aseo, 'tapsmofrq');


if any([isempty(amp_est) isempty(lat_est)])
  doinit = true;
else
  doinit = false;
end

[nsmp_fft, ncomp] = size(erp_fft);            % Number of samples and number of components in frequency domain erp
cfg.aseo.ncomp = ncomp;
cfg.aseo.nsmp_fft = nsmp_fft;
data_init     = real(ifft(data));     % Data in the time domain
data_init     = data_init(1:nsample,:);
data          = data(1:(nsmp_fft/2+1),:); % Get the signal from [0, pi).

%--------------------------------------
% initialization of variables if needed
if isempty(erp_est), erp_est      = erp_fft(1:(nsmp_fft/2+1),:);  end  % ERP waveform, frequency domain
if isempty(amp_est), amp_est      = ones(ntrl, ncomp);    end  % ERP amplitude
if isempty(lat_est), lat_est      = zeros(ntrl, ncomp);   end  % ERP latency
if isempty(noise),   noise        = ones(nsmp_fft/2+1, ntrl); end  %.*repmat(var(data_init,[],1),Nsmp/2+1,1);   % Power spectrum of on-going activity
  
%--------------------------------------------------------------------------------------------
% estimation of the ERP amplitudes and latencies, first pass, starting with uniform estimates 
if doinit,
  amp_in = amp_est;
  lat_in = lat_est;
  for m = 1:ncomp
    [lat_tmp, amp_tmp] = ft_estimate_parameters(cfg, data, erp_est, amp_in, lat_in, m, noise);
    lat_est(:,m) = lat_tmp(:,m);
    amp_est(:,m) = amp_tmp(:,m);
  end
end

%---------------------------------------------------------
% do one iteration of ERP and on-going activity estimation  
for k = 1:numiteration
%reject trials
[rejectflag, corr_est] = ft_rejecttrial(cfg, erp_est, amp_est, lat_est, data);
rejectflag(:) = false;

% ERP estimation in frequency domain
[erp_est, amp_est, lat_est, residual_f] = ft_estimate_erp(cfg, data, erp_est, amp_est, lat_est, noise, rejectflag);

% on-going activity estimation in Time Domain
residual_f     = cat(1, residual_f, conj(residual_f((nsmp_fft/2):-1:2,:)));
residual       = ifft(residual_f,[],1);
residual       = real(residual(1:nsample,:));
residual       = residual - ones(nsample,1)*mean(residual, 1);

temp2 = data_init(1:nsample,:);
temp2 = temp2 - ones(nsample,1)*mean(temp2,1);

temp3 = temp2 - mean(temp2(:,rejectflag==0),2)*ones(1,size(temp2,2));
temp3 = temp3 - ones(size(temp3,1),1)*mean(temp3);
if strcmp(cfg.aseo.noiseEstimate, 'parametric')
  % parametric noise estimate
  
  [noise, ar, sigma] = ft_estimate_ar(residual, rejectflag, Nsmp, maxOrderAR);
  noise              = noise(1:(Nsmp/2+1),:);
  
  % also compute the AR model for the original data
  [orig, ar_orig, sigma_orig] = ft_estimate_ar(temp2, rejectflag, Nsmp, maxOrderAR);
  orig  = orig(1:(Nsmp/2+1),:);
  
  % also compute the AR model for the ERP subtracted data
  [noise3, ar_orig2, sigma_orig2] = ft_estimate_ar(temp3, rejectflag, Nsmp, maxOrderAR);
  noise3 = noise3(1:(Nsmp/2+1),:);
elseif strcmp(cfg.aseo.noiseEstimate, 'non-parametric')
  % non-parametric noise estimate
  
  timeaxis = (1:nsample)./fsample;
  pad      = nsmp_fft./fsample;
  
  [noise, ~, freqoi] = ft_specest_mtmfft(residual',timeaxis,'taper','dpss','pad',pad,'tapsmofrq',tapsmofrq);
  noise              = squeeze(mean(abs(noise).^2))'./(nsample./2);
  ar    = zeros(0,1);
  sigma = [];
  
  orig = ft_specest_mtmfft(temp2',timeaxis,'taper','dpss','pad',pad,'tapsmofrq',tapsmofrq);
  orig = squeeze(mean(abs(orig).^2))'./(nsample./2);
  ar_orig    = zeros(0,1);
  sigma_orig = [];
  
  noise3 = ft_specest_mtmfft(temp3',timeaxis,'taper','dpss','pad',pad,'tapsmofrq',tapsmofrq);
  noise3 = squeeze(mean(abs(noise3).^2))'./(nsample./2);
  ar_orig2    = zeros(0,1);
  sigma_orig2 = [];

else
    error('Please specify cfg.aseo.noiseEstimate as parametric or non-parametric')
end

% adjust the latency such that the mean of latencies is 0.

lat_estMean = mean(lat_est,1);
lat_est     = lat_est - kron(ones(ntrl,1), lat_estMean);

% shift the erp-components along, this assumes to be defined in samples
freqSeq     = 2*pi/nsmp_fft*(0:(nsmp_fft/2))';
for m=1:ncomp
  erp_est(:, m) = erp_est(:,m).*exp(-1i*freqSeq*lat_estMean(m));
end
if ~exist('freqoi','var')
  freqoi = freqSeq;
end

% get time domain version of erp components
tmp_erp_est = cat(1, erp_est, conj(erp_est((nsmp_fft/2):-1:2,:)));
tmp_erp_est = real(ifft(tmp_erp_est));
tmp_erp_est = tmp_erp_est(1:nsample, :);

% scaling
scale       = max(abs(tmp_erp_est));
tmp_erp_est = tmp_erp_est./(ones(size(tmp_erp_est,1),1)*scale);
tmp_amp_est = amp_est.*(ones(ntrl,1)*scale);

output.erp_est    = tmp_erp_est;
output.amp_est    = tmp_amp_est;
output.erp_est_f  = erp_est;
output.lat_est    = lat_est;
output.ar         = ar(:,1);
output.noise      = noise;
output.sigma      = sigma;
output.residual   = residual;
%output.reconstructed = reconstructed;
output.data_init  = data_init(1:nsample,:);
output.rejectflag = rejectflag;
output.corr_est   = corr_est;
output.ar_orig    = ar_orig(:,1);
output.sigma_orig = sigma_orig;
output.ar_orig2    = ar_orig2(:,1);
output.sigma_orig2 = sigma_orig2;
output.rejectflag  = rejectflag;
output.freq_orig   = orig;
[origdata, reconstructed]=ft_reconstruct_erp(output(end));

% figure;
subplot(2,2,1); plot(freqoi,mean(noise,2),freqoi,mean(orig,2),freqoi,mean(noise3,2));xlim([0 60]);drawnow;
subplot(2,2,2); plot([origdata reconstructed]);drawnow;
subplot(2,2,3); plot(tmp_erp_est);drawnow;

[origdata, reconstructed] = ft_reconstruct_erp(output(end),0);
subplot(2,2,4); imagesc(reconstructed'); drawnow;
end
  
params.latency    = lat_est;
params.amplitude  = tmp_amp_est;
params.components = tmp_erp_est;
params.rejectflag = rejectflag;


