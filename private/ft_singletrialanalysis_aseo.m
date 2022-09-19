function [output] = ft_singletrialanalysis_aseo(cfg, data, erp_fft)

% FT_SINGLETRIALANALYSIS_ASEO executes single-trial analysis, using the ASEO
% algorithm (Xu et al, 2009)
%
% Use as
%   [output] = ft_singletrialanalysis_aseo(cfg, data_fft, erp_fft)
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
fsample         = ft_getopt(cfg.aseo, 'fsample');
nsample         = ft_getopt(cfg.aseo, 'nsample'); % the number of original samples in the time domain
amp_est         = ft_getopt(cfg.aseo, 'amp_est', []);
erp_est         = ft_getopt(cfg.aseo, 'erp_est', []);
lat_est         = ft_getopt(cfg.aseo, 'lat_est', []);
noise           = ft_getopt(cfg.aseo, 'noise',   []);
numiteration    = ft_getopt(cfg.aseo, 'numiteration', 1);
tapsmofrq       = ft_getopt(cfg.aseo, 'tapsmofrq');

if all([isempty(amp_est) isempty(lat_est)])
  doinit = true;
else
  doinit = false;
end

[nsmp_fft, ncomp] = size(erp_fft);    % Number of samples and number of components in frequency domain erp
data_init         = real(ifft(data)); % Data in the time domain
data_init         = data_init(1:nsample,:);
data              = data(1:(nsmp_fft/2+1),:); % Get the signal from [0, pi).
ntrl              = size(data,2);

%--------------------------------------
% initialization of variables if needed
if isempty(erp_est), erp_est = erp_fft(1:(nsmp_fft/2+1),:); end  % ERP waveform, frequency domain
if isempty(amp_est), amp_est = ones(ntrl, ncomp);           end  % ERP amplitude
if isempty(lat_est), lat_est = zeros(ntrl, ncomp);          end  % ERP latency
if isempty(noise),   noise   = ones(nsmp_fft/2+1, ntrl);    end  %.*repmat(var(data_init,[],1),nsmp/2+1,1);   % Power spectrum of on-going activity
  
%--------------------------------------------------------------------------------------------
% estimation of the ERP amplitudes and latencies, first pass, starting with uniform estimates
if doinit
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
  %reject trials, based on some heuristics defined in the cfg
  [rejectflag, corr_est] = ft_rejecttrial(cfg, erp_est, amp_est, lat_est, data);
  rejectflag(:)          = false;
  
  % ERP estimation in frequency domain
  [erp_est, amp_est, lat_est, residual_f] = ft_estimate_erp(cfg, data, erp_est, amp_est, lat_est, noise, rejectflag);
  
  % on-going activity estimation in Time Domain
  residual_f = cat(1, residual_f, conj(residual_f((nsmp_fft/2):-1:2,:)));
  residual   = ifft(residual_f,[],1);
  residual   = real(residual(1:nsample,:));
  residual   = residual - ones(nsample,1)*mean(residual, 1);
  
  temp2 = data_init(1:nsample,:);
  temp2 = temp2 - ones(nsample,1)*mean(temp2,1);
  
  temp3 = temp2 - mean(temp2(:,rejectflag==0),2)*ones(1,size(temp2,2));
  temp3 = temp3 - ones(size(temp3,1),1)*mean(temp3);
  if strcmp(cfg.aseo.noiseEstimate, 'parametric')
    % parametric noise estimate
    maxOrderAR = ft_getopt(cfg.aseo, 'maxOrderAR', 5);

    [noise, ar, sigma] = ft_estimate_ar(residual, rejectflag, nsmp_fft, maxOrderAR);
    noise              = noise(1:(nsmp_fft/2+1),:);
    
    % also compute the AR model for the original data
    [orig, ar_orig, sigma_orig] = ft_estimate_ar(temp2, rejectflag, nsmp_fft, maxOrderAR);
    orig                        = orig(1:(nsmp_fft/2+1),:);
    
    % also compute the AR model for the ERP subtracted data
    [noise2, ar_orig2, sigma_orig2] = ft_estimate_ar(temp3, rejectflag, nsmp_fft, maxOrderAR);
    noise2                          = noise2(1:(nsmp_fft/2+1),:);
    
  elseif strcmp(cfg.aseo.noiseEstimate, 'nonparametric')
    % non-parametric noise estimate
    timeaxis = (1:nsample)./fsample;
    pad      = nsmp_fft./fsample;
    optarg   = {'taper','dpss','pad',pad,'tapsmofrq',tapsmofrq};
    
    [noise, dum, freqoi] = ft_specest_mtmfft(residual',timeaxis,optarg{:});
    noise                = squeeze(mean(abs(noise).^2))'./(nsample./2);
    
    orig = ft_specest_mtmfft(temp2',timeaxis,optarg{:});
    orig = squeeze(mean(abs(orig).^2))'./(nsample./2);
    
    noise2 = ft_specest_mtmfft(temp3',timeaxis,optarg{:});
    noise2 = squeeze(mean(abs(noise2).^2))'./(nsample./2);
    
  else
    error('Please specify cfg.aseo.noiseEstimate as parametric or nonparametric')
  end
  
  % adjust the latency such that the mean of latencies is 0.
  lat_estMean = mean(lat_est,1);
  lat_est     = lat_est - kron(ones(ntrl,1), lat_estMean);
  
  % shift the erp-components along, this assumes to be defined in samples
  freqSeq = 2*pi/nsmp_fft*(0:(nsmp_fft/2))';
  for m = 1:ncomp
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
  
  output(k).erp_est     = tmp_erp_est;
  output(k).amp_est     = tmp_amp_est;
  output(k).erp_est_f   = erp_est;
  output(k).lat_est     = lat_est;
  output(k).noise       = noise;
  output(k).residual    = residual;
  output(k).data_init   = data_init(1:nsample,:);
  output(k).rejectflag  = rejectflag;
  output(k).corr_est    = corr_est;
  output(k).rejectflag  = rejectflag;
  output(k).freq_orig   = orig;
  if exist('ar', 'var')
    output(k).ar = ar(:,1);
    output(k).sigma       = sigma;
    output(k).ar_orig     = ar_orig(:,1);
    output(k).sigma_orig  = sigma_orig;
    output(k).ar_orig2    = ar_orig2(:,1);
    output(k).sigma_orig2 = sigma_orig2;
  end
  [origdata_avg, reconstructed_avg] = ft_reconstruct_erp(output(end));
  
  figure;
  subplot(2,2,1); plot(freqoi,mean(noise,2),freqoi,mean(orig,2),freqoi,mean(noise2,2));xlim([0 60]);drawnow;
  subplot(2,2,2); plot([origdata_avg reconstructed_avg]);drawnow;
  subplot(2,2,3); plot(tmp_erp_est);drawnow;
  
  [origdata, reconstructed] = ft_reconstruct_erp(output(k),0);
  subplot(2,2,4); imagesc(reconstructed'); drawnow;
end

function  [power, coeffAR, sigma] = ft_estimate_ar(dat, rejectflag, nfft, maxOrderAR)

% Use Levinson-Durbin to estimate the AR coefficicent
% Use Bayesian Information Criterion for AR order determination
%
% --Input
% dat        : On-going activity in time-domain
% rejectflag : Indicating the trials that will be rejected. If would not
%              reject any trials, set rejectflag=[0 0 0 0 ....] .
% nfft       : Sample number to represent the power spectrum
%
%--Output
% power   : Power spectrum of the estimated AR model
% coeffAR : Estimated AR coefficients
% sigma   : Power of the input white noise in the AR model

[dum, nrpt_orig] = size(dat);

% reject the poor trials
dat          = dat(:,~rejectflag);
[nsmp, nrpt] = size(dat);

coeffARTrial = zeros(maxOrderAR, maxOrderAR);
costFun      = zeros(maxOrderAR,1);

% calculate the autocorrelation
ruo0 = sum(sum(dat.*conj(dat), 1), 2)/nrpt/nsmp;
ruo  = zeros(maxOrderAR,1);
for k=1:maxOrderAR
    ruo(k) = sum(sum(dat(1:end-k,:) .* conj(dat(k+1:end,:)), 1),2)/nrpt/(nsmp-k);
end

% initialize
kai               = -ruo(1)/ruo0;
theta             = kai;
coeffARTrial(1,1) = kai;
sigmaSquare       = ruo0-abs(ruo(1)).^2/ruo0;

sigma             = zeros(1, maxOrderAR);
sigma(1)          = sigmaSquare;
costFun(1)        = nrpt*nsmp.*log(sigmaSquare)+4;
for n = 1:maxOrderAR-1
  kai         = -1*(ruo(n+1) + ruo(n:-1:1)'* theta)/sigmaSquare;
  sigmaSquare = sigmaSquare*(1-abs(kai).^2);
  theta       = [theta; 0] + kai*[theta(n:-1:1); 1 ];
  
  coeffARTrial(1:n+1,n+1) = theta;
  sigma(n+1)              = sigmaSquare;
  costFun(n+1)            = nrpt*nsmp.*log(sigmaSquare)+(n)*log(nrpt*nsmp); % BIC cost function
end

[dum, orderAR] = min(real(costFun));
fprintf('the optimal AR-model order is %d\n', orderAR);
coeffAR      = coeffARTrial(1:orderAR, orderAR);

% calculate the frequency response of the AR model
temp    = ones(nfft, 1);
freqSeq = 2*pi/nfft*(0:(nfft-1))';
for k = 1:orderAR
  temp = temp + coeffAR(k).*exp(-1i.*freqSeq.*k);
end

% calculate the input noise power
X = dat(orderAR+1:end,:);
for k = 1:orderAR
    X = X + coeffAR(k)*dat(orderAR+1-k:end-k,:);
end
sigma = sum(X.*conj(X),1)./(nsmp-orderAR);

% calculate the power spectrum
power   = (1./ (abs(temp).^2))*sigma;
coeffAR = -1*kron(coeffAR, ones(1, nrpt));

% create output
temp      = mean(power,2);
powertemp = temp*ones(1,nrpt_orig);
powertemp(:,~rejectflag) = power;
power     = powertemp;

temp      = mean(sigma,2);
sigmatemp = temp*ones(1,nrpt_orig);
sigmatemp(:,~rejectflag) = sigma;
sigma     = sigmatemp;

function [orig, reconstructed] = ft_reconstruct_erp(input, avgflag)

if nargin<2
  avgflag = true;
end

[nsmp,ncomp]  = size(input.erp_est);
ntrl          = size(input.amp_est,1);
rejectflag    = input.rejectflag;

if avgflag
  reconstructed = zeros(nsmp,1);
  orig          = zeros(nsmp,1);
  for k = 1:ntrl
    if ~rejectflag(k)
      orig = orig + input.data_init(:,k);
      
      for m=1:ncomp
        tmp           = input.erp_est(:,m);
        reconstructed = reconstructed + input.amp_est(k,m).*fun_shift(tmp, input.lat_est(k,m), 1);
      end
    end
  end
  reconstructed = reconstructed./sum(rejectflag==0);
  orig          = orig./sum(rejectflag==0);
else
  reconstructed = zeros(nsmp,sum(rejectflag==0));
  orig          = input.data_init(:,rejectflag==0);
  indx = 0;
  for k = 1:ntrl
    if ~rejectflag(k)
      indx = indx+1;
      for m=1:ncomp
        tmp                   = input.erp_est(:,m);
        reconstructed(:,indx) = reconstructed(:,indx) + input.amp_est(k,m).*fun_shift(tmp, input.lat_est(k,m), 1);
      end
    end
  end
end

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
fsample         = ft_getopt(cfg.aseo, 'fsample');
jitter          = ft_getopt(cfg.aseo, 'jitter'); % Latency search window defined in s

inv_noise           = 1./noise;
data([1 end],:)     = data([1 end],:)/2;     % JM note: this probably has something to do with DC and Nyquist
erp_est([1 end], :) = erp_est([1 end], :)/2; % JM note: this probably has something to do with DC and Nyquist
[nsmp, Ncomp]       = size(erp_est);
ntrl                = size(data,2);
freqSeq             = 2*pi/(2*nsmp-2)*(0:(nsmp-1))';
  
% Estimate ERP latencies and amplitudes trial by trial
sel       = [1:(compNo-1) (compNo+1):size(amp_est,2)];%setdiff(1:size(amp_est,2), compNo);
fft_point = round(2*nsmp-2);
index1 = round(-jitter(compNo).*fsample + fft_point/2);
index2 = round( jitter(compNo).*fsample + fft_point/2);
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
zero_c(  1:index1,:) = false;
zero_c(index2:end,:) = false;
[i1,i2] = find(zero_c);
i1      = i1 - round(fft_point/2);

cnt = 0;
for trialNo = 1:ntrl
  tmp_i1 = i1(i2==trialNo);
  if isempty(tmp_i1)
    cnt = cnt+1;
    % fall back to the original, when no local maximum is found in the
    % search window, this then probably takes an edge
    ruo        = T(index1:index2,trialNo);
    index      = round(numel(ruo)./2);
    lat_est(trialNo, compNo) = (-jitter(compNo) + (index-1));%*searchGrid)*(fsample/1000); %searchWindow(index) in samples (as opposed as original implementation;%
  else
    [dum, index] = min(abs(tmp_i1));
    index      = tmp_i1(index);
    lat_est(trialNo, compNo) = (index-1); %searchWindow(index) in samples (as opposed as original implementation;
  end
end
fprintf('the number of trials for which no max was found = %d\n',cnt);

% Estimate the ERP amplitudes using Least-Squares in the frequency domain
for trialNo = 1:ntrl
  A     = exp(-1i*freqSeq* lat_est(trialNo,:) ).*erp_est;
  A_tmp = A.*(inv_noise(:,trialNo) * ones(1, Ncomp));
  x     = data(:,trialNo);
  amp_estTmp = real(A'*A_tmp)\real(A_tmp'*x);
  
  amp_estTmp(amp_estTmp<0) = nan; % disallow for polarity flips
  amp_est(trialNo, :)      = amp_estTmp.';
end

if any(~isfinite(amp_est(:)))
  for k = 1:size(amp_est,2)
    amp_est(~isfinite(amp_est(:,k)),k) = nanmedian(amp_est(:,k));%./10;
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

function dT = ft_deriv(T)

% derivative along the columns, where the 2:end-1 elements are the average
% of the n-1 and n+1 differences

[i1, i2]  = size(T);
dT(i1,i2) = 0;
dT(1,  :) = T(2,:)  - T(1,:);
dT(i1, :) = T(i1,:) - T(i1-1,:);
dT(2:(i1-1),:) = (T(3:i1,:)-T(1:(i1-2),:))./2;
  
function [seq] = fun_shift(seq, step, dim)

% step < 0    -----  Left or up
% step > 0    -----  Right or Down
% Dim =1      -----  Row
% Dim =2      -----  Col
[row, col] = size(seq);
step       = -1*round(step);
    
if (dim==2) && (abs(step)>=col)
  seq = zeros(row, col);
  return
end
    
if (dim~=2) && (abs(step) >= row)
  seq = zeros(row, col);
  return
end
    
if step>0
  if dim ==2
    seq= [seq(:, step+1:end) zeros(row, step)];
  else
    seq= [seq(step+1:end,:); zeros(step, col)];
  end
end
    
if step <0
  if dim ==2
    seq= [zeros(row, abs(step)) seq(:, 1:end+step) ];
  else
    seq= [zeros(abs(step), col); seq(1:end+step,:) ];
  end
end

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
  
ntrl          = size(data,2);
ncomp         = size(amp_in,2);
nsmp          = size(erp_in, 1); % Sample number and component number
freqSeq       = 2*pi/(2*nsmp-2)*(0:(nsmp-1))';
  
invsqrt_noise = 1./sqrt(noise);
acceptIndex   = find(~rejectflag);
erp_est       = zeros(nsmp,ncomp);
for k = 1:nsmp
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
  lat_est(:,compNo) = tmp_lat(:,compNo);
  amp_est(:,compNo) = tmp_amp(:,compNo);
end
      
% Compute residual signal by removing ERPs
residual_f      = zeros(nsmp, ntrl);
reconstructed_f = zeros(nsmp, ntrl);
for trialNo = 1:ntrl
  reconstructed_f(:,trialNo) = ( exp(-1i*freqSeq*lat_est(trialNo,:) ).*erp_est )* amp_est(trialNo, :).';
  residual_f(:,trialNo)      = data(:,trialNo) - reconstructed_f(:,trialNo);
end

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
ncomp   = size(amp_est,2);
ntrl    = size(data,2);


nsmp = size(data,1);
freqSeq      = 2*pi/(2*nsmp-1)*(0:(nsmp-1))';
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
