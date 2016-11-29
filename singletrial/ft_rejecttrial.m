function [rejectflag, corr_est] = ft_rejecttrial(cfg, erp_est, amp_est, lat_est, data)

% Description: Reject trials
%
% INPUT----
% erp_est: Estimates of ERP waveforms in frequency domain
% amp_est: Estimates of ERP components
% lat_est: Estimates of ERP components
% data:    Observed data samples in frequency domain
%
% OUTPUT----
% rejectflag : Boolean vector indicating whether or not a trial is rejected
% corr_est   : Correlation between the original data and the recovered signal

thresholdAmpH   = ft_getopt(cfg.aseo, 'thresholdAmpH'); % maximum acceptable amplitude of trials, times of avergae amplitude
thresholdAmpL   = ft_getopt(cfg.aseo, 'thresholdAmpL'); % minimum  acceptable amplitude of trials, times of avergae amplitude
thresholdCorr   = ft_getopt(cfg.aseo, 'thresholdCorr'); % minimum correlation with the original data
ncomp   = ft_getopt(cfg.aseo, 'ncomp'); 
ntrl    = ft_getopt(cfg.aseo, 'ntrl');


Nsmp = size(data,1); 
freqSeq      = 2*pi/(2*Nsmp-1)*(0:(Nsmp-1))';
rejectflag   = false(ntrl,1);
corr_est     = zeros(ntrl,1);

for k = 1:ncomp
  rejectflag = rejectflag | sign(amp_est(:,k))~=median(sign(amp_est(:,k)));
end

% for k = 1:ncomp
%   amp_est_k    = amp_est(:,k);
%   sel          = 1:numel(amp_est_k);%sign(amp_est_k)==polarity(k);
%   amp_est_mean = mean(amp_est_k(sel));
%   
%   % check the amplitudes for each trial, but take the polarity into account
%   if ~isempty(polarity)
%     switch polarity(k)
%       case 1
%         rejectflag = rejectflag | (amp_est_k>thresholdAmpH*amp_est_mean) | (amp_est_k<thresholdAmpL*amp_est_mean); 
%       case -1
%         rejectflag = rejectflag | (amp_est_k<thresholdAmpH*amp_est_mean) | (amp_est_k>thresholdAmpL*amp_est_mean);
%     end
%   else
%   
%     rejectflag = rejectflag|(amp_est_k<thresholdAmpL*amp_est_mean)|(amp_est_k>thresholdAmpH*amp_est_mean);
%   end
% end

% check the correlation between the original data samples and estimated ERPs
for k = 1:ntrl
  signal = exp(-1i*freqSeq* lat_est(k,:) ).*erp_est;     
  signal = signal*amp_est(k, :).';
        
  corr_tmp      = real(signal'*data(:,k))/norm(signal)/norm(data(:,k)) ; 
  rejectflag(k) = rejectflag(k) | real(corr_tmp) < thresholdCorr;
  corr_est(k)   = corr_tmp;
end
