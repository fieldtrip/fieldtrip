function test_ft_crossfrequencyanalysis

% MEM 2gb
% WALLTIME 00:15:00

% TEST ft_crossfrequencyanalysis

clear all
close all

%% generate simulation data

%  channels
N  = 2;
s = zeros(N,4000);
for  i = 1:N
  num     = 45;                   % number of alpha cycles
  fs      = 1000;                 % sampling frequency
  hf      = 70;                   % gamma frequency
  shift   ='lead'; timdiff = 10;  % for directionality only
  a = 10; c = 6;                  % a(slope)/c(threshold)  :   Sigmoid function parameter sigmf(x, [a, c]) = 1./(1 + exp(-a*(x-c)))
  n1      = rand(1);              % Gaussian white noise level
  n2      = rand(1);              % pink noise level
  [sig,T] = inhibition(num,fs,shift,timdiff,hf,a,c,n1,n2);
  s(i,:)  = sig(4,1:4000);
end

M = 20;
ftdata = zeros(M,N,4000);  %  M=trials * N=Channels * times

for  j = 1:M
  ftdata(j,:,:) = s+rand(N,4000);
end

% represent it in FieldTrip fashion
data              = [];
data.time         = cell(1,M);
data.trial        = cell(1,M);
for i =1:M
  data.time{1,i}   = T(1:4000);
  data.trial{1,i}  = squeeze(ftdata(i,:,:));
end
data.fsample       = 1000;
data.label         = cell(N,1);

for  j = 1:N
  data.label{j,1} = strcat('chan',num2str(j));
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = 4:1:20;       % interest low frequency range of CFC
f2 = 30:10:150;    % interest high frequency range of CFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract low frquency signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg              = [];
cfg.output       = 'fourier';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          =  f1;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;
cfg.toi          = 0.5:1/data.fsample:3.5;
LFsig            = ft_freqanalysis(cfg, data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract high frequency evelope signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg              = [];
cfg.output       = 'fourier';
cfg.channel      = 'all';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          =  f2;
cfg.t_ftimwin    = 5./cfg.foi;
cfg.toi          = 0.5:1/data.fsample:3.5;
HFsig            = ft_freqanalysis(cfg, data);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the actual testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cfg              = [];
cfg.method       = 'coh';
cfg.keeptrials   = 'no';
CFC              = ft_crossfrequencyanalysis(cfg, LFsig, HFsig);

subplot(411)
MI = squeeze(CFC.crsspctrm(1,:,:));
imagesc(f1, f2, MI');
% set(gca, 'Fontsize',20)
axis xy;
xlabel('Low frequency  (Hz)');
ylabel('High frequency (Hz)');
title('Coherence')
axis xy; colorbar

cfg              = [];
cfg.method       = 'plv';
cfg.keeptrials   = 'no';
CFC              = ft_crossfrequencyanalysis(cfg, LFsig, HFsig);

subplot(412)
MI = squeeze(CFC.crsspctrm(1,:,:));
imagesc(f1, f2, MI');
% set(gca, 'Fontsize',20)
axis xy;
xlabel('Low frequency  (Hz)');
ylabel('High frequency (Hz)');
title('Phase locking value')
axis xy; colorbar

cfg              = [];
cfg.method       = 'mvl';
cfg.keeptrials   = 'no';
CFC              = ft_crossfrequencyanalysis(cfg,LFsig,HFsig);

subplot(413)
MI = squeeze(CFC.crsspctrm(1,:,:));
imagesc(f1, f2, MI');
% set(gca, 'Fontsize',20)
axis xy;
xlabel('Low frequency  (Hz)');
ylabel('High frequency (Hz)');
title('mean vector length')
axis xy; colorbar

cfg              = [];
cfg.method       = 'mi';
cfg.keeptrials   = 'no';
CFC              = ft_crossfrequencyanalysis(cfg,LFsig,HFsig);

subplot(414)
MI = squeeze(CFC.crsspctrm(1,:,:));
imagesc(f1, f2, MI');
% set(gca, 'Fontsize',20)
axis xy;
xlabel('Low frequency  (Hz)');
ylabel('High frequency (Hz)');
title('Modulation index')
axis xy; colorbar

end % main function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigs, T2] = inhibition(num,fs,shift,timdiff,hf,a,c,n1,n2)

% generate simulation data
% input:
% num:     number of alpha cycles
% fs:      sampling frequency
% hf:      gamma frequency
% shift:   create directionality either alpha leads gamma or alpha lags gamma
% timdiff:  alpha and gamma time lag
% a(slope)/c(threshold)  :   Sigmoid function parameter sigmf(x, [a, c]) = 1./(1 + exp(-a*(x-c)))
%  n1                    :   Gaussian white noise level
%  n2                    :   pinck noise level
%
% output:
% sigs:     four signals {alpha0;alphas;gamma;mix}  4*T2   we are looking CFC at mix channel sigs(4,:)
% T2        time series index

% copyright @Haiteng Jiang

dt         = 1/fs;
Time       =  0.1 + 0.02.*randn(num,1);  % alpha range cycle length 80-120 ms
amps       =  2 + 0.5.*randn(num,1);     % fluctuated alpha amplitude
alpha      = [];
alpha2     = [];

% generate fluctuated amplitude alpha signal
for i = 1:num
  T    = Time(i);
  t    = 0:dt:T-dt;
  N    = length(t);
  sig  = (1+sin(1/T* 2*pi*t+1.5*pi));
  alpha2 = [alpha2 sig];
  sig   = amps(i)*sig;       % fluctuated alpha amplitude
  alpha = [alpha sig];
end

alpha      = 6* alpha;     % enhanced amplitude to be more real
T1         = 0:dt:length(alpha)*dt-dt;
T2         = 0:dt:(length(alpha)-timdiff)*dt-dt;
%  sigmoid threshold gamma
gammas     = (1-1./(1 + exp(-a*(alpha-c)))).*(sin(2*pi*hf*T1)+1);

% shift to creat directionality not mean for CFC
switch shift
  case {'lead'}       % alpha leads gamma
    alphastemp = zeros(1,length(alpha));
    for i = 1:length(alpha)-timdiff
      alphastemp(i) = alpha(i+timdiff);
    end
    alphas = alphastemp(1:end-timdiff);
    gamma  = gammas(1:end-timdiff);
    
  case{'delay'}      % alpha  lags gamma
    alphastemp = zeros(1,length(alpha));
    for i = timdiff+1:length(alpha)
      alphastemp(i) = alpha(i-timdiff);
    end
    alphas(1:length(T2))= alphastemp(timdiff+1:length(alpha));
    gamma  = gammas(timdiff+1:end);
  otherwise         % no directionality
    alphas = alpha(1:end-timdiff);
    gamma  = gammas(1:end-timdiff);
end

alpha0  = alpha2(1:length(T2));
mix    =  alphas+gamma+n1*randn([1,length(T2)])+n2*pinknoise(numel(T2)); % signal we analysis
sigs = [alpha0;alphas;gamma;mix];

  function x = pinknoise(Nx)
    % pink noise
    B = [0.049922035 -0.095993537 0.050612699 -0.004408786];
    A = [1 -2.494956002   2.017265875  -0.522189400];
    nT60 = round(log(1000)/(1-max(abs(roots(A))))); % T60 est.
    v = randn(1,Nx+nT60); % Gaussian white noise: N(0,1)
    x = filter(B,A,v);    % Apply 1/F roll-off to PSD
    x = x(nT60+1:end);    % Skip transient response
  end

end
