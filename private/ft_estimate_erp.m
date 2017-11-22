function [erp_est, amp_est, lat_est, residual_f] = ft_estimate_erp(cfg, data, erp_in, amp_in, lat_in, noise, rejectflag)

% Estimate the ERP waveforms, amplitudes and latencies in the frequency domain
%
% INPUT----
% data:  Observed data samples in frequency domain
% erp_est    : The most-recent estimates of the ERP waveforms in frequency domain
% amp_est    : The most-recent estimates of ERP amplitudes
% lat_est    : The most-recent estimates of ERP latencies
% noise      : The most=recent estimates of spectrum of on-going activity
% rejectflag : Reject flags
%
% OUTPUT ----
% erp_est    : Updated estimates of the ERP waveforms in frequency domain
% amp_est    : Updated estimates of ERP amplitudes
% lat_est    : Updated estimates of ERP latencies
% residual_f : Residual signal after removing ERPs in frequency domain
  
ntrl          = ft_getopt(cfg.aseo, 'ntrl');
ncomp         = ft_getopt(cfg.aseo, 'ncomp');
Nsmp          = size(erp_in, 1); % Sample number and component number
freqSeq       = 2*pi/(2*Nsmp-2)*(0:(Nsmp-1))';
  
invsqrt_noise = 1./sqrt(noise);
acceptIndex   = find(~rejectflag);
for k = 1:Nsmp
  A_tilde = amp_in(acceptIndex,:).*exp(-1i*freqSeq(k)*lat_in(acceptIndex,:)).* (invsqrt_noise(k, acceptIndex).'*ones(1,ncomp)) ;  
  X_tilde = data(k,acceptIndex).'.*invsqrt_noise(k,acceptIndex).';
  S_tilde = (A_tilde'* A_tilde)\(A_tilde'*X_tilde);
       
  erp_est(k,:) = S_tilde.';   
end
     
% update the latency and amplitude estimates
lat_est = lat_in;
amp_est = amp_in;
for compNo = 1:ncomp
  [tmp_lat, tmp_amp] = ft_estimate_parameters(cfg, data, erp_est, amp_in, lat_in, compNo, noise);
  lat_est(:,compNo) = tmp_lat(:,compNo);%+lat_in(:,compNo);
  amp_est(:,compNo) = tmp_amp(:,compNo);
end
      
% Compute residual signal by removing ERPs
residual_f      = zeros(Nsmp, ntrl);
reconstructed_f = zeros(Nsmp, ntrl); 
for trialNo = 1:ntrl
  reconstructed_f(:,trialNo) = ( exp(-1i*freqSeq*lat_est(trialNo,:) ).*erp_est )* amp_est(trialNo, :).';
  residual_f(:,trialNo) = data(:,trialNo) - reconstructed_f(:,trialNo);
end
