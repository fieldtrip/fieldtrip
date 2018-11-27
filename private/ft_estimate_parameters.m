function [lat_est, amp_est] = ft_estimate_parameters(cfg, data, erp_est, amp_est, lat_est, compNo, noise)

% Estimate the ERP amplitudes and latencies in frequency domain
%
% INPUT----
% data    : Observed data samples in frequency domain
% erp_est : The most-recent estimates of the ERP waveforms in frequency domain
% amp_est : The most-recent estimates of ERP amplitudes
% lat_est : The most-recent estimates of ERP latencies
% compNo  : Component # 
% noise   : The most=recent estimates of spectrum of on-going activity
%
% OUTPUT ----
% amp_est     : Updated estimates of ERP amplitudes
% lat_est     : Updated estimates of ERP latencies

% define some options locally
period          = ft_getopt(cfg.aseo, 'sampPeri');
jitter = ft_getopt(cfg.aseo, 'jitter'); % Latency search window defined in main_ASEO.m
searchGrid      = ft_getopt(cfg.aseo, 'searchGrid'); % Latency search step
ntrl            = ft_getopt(cfg.aseo, 'ntrl');
fsample         = ft_getopt(cfg.aseo, 'fsample');
%{
maxOrderAR      = ft_getopt(cfg.aseo, 'maxOrderAR', 10);
thresholdAmpH   = ft_getopt(cfg.aseo, 'thresholdAmpH'); % maximum acceptable amplitude of trials, times of avergae amplitude
thresholdAmpL   = ft_getopt(cfg.aseo, 'thresholdAmpL'); % minimum  acceptable amplitude of trials, times of avergae amplitude
thresholdCorr   = ft_getopt(cfg.aseo, 'thresholdCorr'); % minimum correlation with the original data
nsample         = ft_getopt(cfg.aseo, 'nsample'); % the number of original samples in the time domain
amp_est         = ft_getopt(cfg.aseo, 'amp_est', []);
erp_est         = ft_getopt(cfg.aseo, 'erp_est', []);
lat_est         = ft_getopt(cfg.aseo, 'lat_est', []);
noise           = ft_getopt(cfg.aseo, 'noise', []);
numiteration    = ft_getopt(cfg.aseo, 'numiteration', 1);
ntrl            = ft_getopt(cfg.aseo, 'ntrl', size(data_fft, 2));
%}


inv_noise   = 1./noise;
  
data([1 end],:) = data([1 end],:)/2;     % JM note: this probably has something to do with DC and Nyquist
erp_est([1 end], :) = erp_est([1 end], :)/2; % JM note: this probably has something to do with DC and Nyquist
[Nsmp, Ncomp]       = size(erp_est);
freqSeq             = 2*pi/(2*Nsmp-2)*(0:(Nsmp-1))';
  
% Estimate ERP latencies and amplitudes trial by trial 
sel       = [1:(compNo-1) (compNo+1):size(amp_est,2)];%setdiff(1:size(amp_est,2), compNo);
fft_point = round((2*Nsmp-2)*period/searchGrid);
searchWindow_index1 = round(jitter(compNo,1)/searchGrid) + round(fft_point/2);
searchWindow_index2 = round(jitter(compNo,2)/searchGrid) + round(fft_point/2);    
temp = zeros(size(noise));
for trialNo = 1:ntrl
             
  amp_estTmp = amp_est(trialNo, sel);
  lat_estTmp = lat_est(trialNo, sel);
  erp_estTmp = erp_est(:,       sel);
  
  % Signal after removal of the other ERP components
  signal = data(:,trialNo) - (exp(-1i*freqSeq*lat_estTmp ).*erp_estTmp) * amp_estTmp.';  
           
  % Estimate the ERP latency
  temp(:,trialNo) = inv_noise(:,trialNo).*conj(signal).*(erp_est(:,compNo));
end  
temp_tau  = fft(temp, fft_point);  % Do Fourier transform on all trials at once
temp_tau  = fftshift(temp_tau, 1); % and shift

% NOTE: the IEEE paper mentions an ifft on
% inv_noise.*signal.*conj(erp_est), however
% fftshift(fft(x)) = ifftshift(ifft(conj(x))), so this is equivalent

% now extract the latencies in a different way than the original
% implementation
T  = real(temp_tau);
dT = ft_deriv(T); 

% zero crossings of the downward sloping derivative
zero_c = sign(dT(1:end-1,:))>sign(dT(2:end,:));
zero_c(1:(round(fft_point/2) + jitter(compNo,1)),:)     = false;
zero_c(  (round(fft_point/2) + jitter(compNo,2)):end,:) = false;
[i1,i2] = find(zero_c);
i1      = i1 - round(fft_point/2);
cnt = 0;
for trialNo = 1:ntrl
  tmp_i1 = i1(i2==trialNo);
  if isempty(tmp_i1)
    cnt = cnt+1;
    % fall back to the original, when no local maximum is found in the
    % search window, this then probably takes an edge
    ruo        = real(temp_tau(searchWindow_index1:searchWindow_index2,trialNo));
    [~, index] = max(ruo,[],1);
    index = round(numel(ruo)./2);
    lat_est(trialNo, compNo) = (jitter(compNo,1) + (index-1)*searchGrid)*(fsample/1000); %searchWindow(index) in samples (as opposed as original implementation;%  
  else
    [~, index] = min(abs(tmp_i1));
    index      = tmp_i1(index);
    lat_est(trialNo, compNo) = ((index-1)*searchGrid)*(fsample/1000); %searchWindow(index) in samples (as opposed as original implementation;  
  end
end
fprintf('the number of trials for which no max was found = %d\n',cnt);

% Estimate the ERP amplitudes using Least-Squares in the frequency domain
for trialNo = 1:ntrl
  A     = exp(-1i*freqSeq* lat_est(trialNo,:) ).*erp_est; 
  A_tmp = A.*(inv_noise(:,trialNo) * ones(1, Ncomp)) ;
  x     = data(:,trialNo);
  amp_estTmp          = real(A'*A_tmp)\real(A_tmp'*x);
  
  %amp_estTmp = real(A(:,compNo)'*A_tmp(:,compNo))\real(A_tmp(:,compNo)'*x);
  
  amp_estTmp(amp_estTmp<0) = nan;
  
  amp_est(trialNo, :) = amp_estTmp.';
  
  %amp_est(trialNo, compNo) = amp_estTmp.'; % should the amplitude of all components be updated?
end

if any(~isfinite(amp_est(:)))
%     keyboard
for k = 1:size(amp_est,2)
  amp_est(~isfinite(amp_est(:,k)),k) = nanmedian(amp_est(:,k))./10;
end
end

% erp_tmp = real(ifft(cat(1,erp_est,conj(erp_est(end-1:-1:2,:)))));
% dat_tmp = real(ifft(cat(1,data,conj(data(end-1:-1:2,:)))));  
% 
% erp_tmp = erp_tmp(1:options.nsample,:);
% dat_tmp = dat_tmp(1:options.nsample,:);
% dat_tmp_new = zeros(size(dat_tmp));
% 
% for trialNo = 1:Nrpt
%   % shift dat_tmp in the time domain to match the latency, negative sign is
%   % on purpose
%   dat_tmp_new(:,trialNo) = fun_shift(dat_tmp(:,trialNo), -lat_est(trialNo, compNo), 1);
% end
%   
% % do a regression
% X   = erp_tmp(:,compNo)-mean(erp_tmp(:,compNo));
% X   = X./max(X);
% W   = diag((X-median(X)).^2);
% tmp = (X'*W*X)\X'*W*dat_tmp_new;
% 
% %tmp = dat_tmp_new'/[X ones(options.nsample,1)]';
% 
% amp_est(:,compNo) = tmp';
% 

function dT = ft_deriv(T)

% derivative along the columns, where the 2:end-1 elements are the average
% of the n-1 and n+1 differences

[i1,i2]   = size(T);
dT(i1,i2) = 0;

dT(1,:)  = T(2,:)  - T(1,:);
dT(i1,:) = T(i1,:) - T(i1-1,:);
dT(2:(i1-1),:) = (T(3:i1,:)-T(1:(i1-2),:))./2;
  
