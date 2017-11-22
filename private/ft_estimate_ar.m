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

[~, nrpt_orig] = size(dat);

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

[~, orderAR] = min(real(costFun));
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
end;
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
