function test_di

% MEM 2gb
% WALLTIME 00:20:00
% DEPENDENCY ft_connectivityanalysis ft_connectivity_mutualinformation 

%% Simulate two situations as in analogy to Robin's paper

% situation 1: di source->target and dfi 
% situation 2: di source->target, no dfi
nstim    = 500; 
isi      = ceil(200.*rand(1,nstim));
amp_smp  = cumsum(isi);
amp_stim = rand(1,nstim);

feature1         = zeros(1, max(amp_smp+150));
feature1(amp_smp) = amp_stim;

feature2a = zeros(1, max(amp_smp+150));
feature2b = zeros(1, max(amp_smp+150));
feature2a(amp_smp(amp_stim<0.5)) = amp_stim(amp_stim<0.5).*2;
feature2b(amp_smp(amp_stim>0.5)) = amp_stim(amp_stim>0.5);

krn1   = cat(1,zeros(50,1),gausswin(50,5)); krn1 = krn1./sum(krn1);
krn2   = cat(1,zeros(20,1),gausswin(50,5)); krn2 = krn2./sum(krn2);
krn2b  = cat(1,zeros(70,1),gausswin(50,5)); krn2b = krn2b./sum(krn2b);
krn2c  = gausswin(50,5);krn2c = krn2c./sum(krn2c);

% situation 1
source1 = conv(feature1, krn1', 'same') + randn(1,numel(feature1))./100;
target1 = conv(source1, krn2', 'same') + randn(1,numel(feature1))./100;

% situation 2
nstim2    = 250; 
isi2      = ceil(400.*rand(1,nstim2));
amp_smp2  = cumsum(isi2); amp_smp2(amp_smp2>amp_smp(end)+150) = [];
nstim2    = numel(amp_smp2);
amp_stim2 = rand(1,nstim2); 

noisex  = zeros(1, numel(feature1));
noisex(amp_smp2) = amp_stim2;

noise2t  = zeros(size(noisex));
noise2s  = zeros(size(noisex));
for k = 1:numel(amp_smp)
  noise2t = conv(noisex, krn2c', 'same') +randn(1,numel(noisex))./250;
  noise2s = conv(noisex, krn2c', 'same') +randn(1,numel(noisex))./250;
end
  

source2 = conv(feature2a, krn1', 'same') ; %+ conv(noise2s,krn2c','same')  + randn(1,numel(feature2b))./200;
target2 = conv(feature2b, krn2b', 'same'); %+ conv(noise2t,krn2c','same') + randn(1,numel(feature2b))./200;
source2(1:end-10)   = source2(1:end-10)+noise2s(1:end-10);
target2(11:end) = target2(11:end)+noise2t(1:end-10);

fs    = 100;        % in Hz
tim   = (1:numel(feature1))./fs;

data1          = [];
data1.trial{1} = [feature1; source1; target1];
data1.time{1}  = tim;
data1.label    = {'feature';'source';'target'};
data1.fsample  = fs;

data2 = data1;
data2.trial{1} = [feature2a+feature2b; source2; target2];

%% Run
cfg         = [];
cfg.method  = 'mi';
cfg.refindx = 'all';
cfg.mi.lags = (-fs/2:2:fs/2)./fs;
%mi1 = ft_connectivityanalysis(cfg, data1);
%mi2 = ft_connectivityanalysis(cfg, data2);
cfg.mi.precondition = true;
mi1b = ft_connectivityanalysis(cfg, data1);
mi2b = ft_connectivityanalysis(cfg, data2);

cfg            = [];
cfg.method     = 'di';
cfg.refindx    = 'all';
cfg.di.lags    = (1:1:fs/2)./fs;
%di1 = ft_connectivityanalysis(cfg, data1);
%di2 = ft_connectivityanalysis(cfg, data2);
cfg.di.precondition = true;
di1b = ft_connectivityanalysis(cfg, data1);
di2b = ft_connectivityanalysis(cfg, data2);

cfg = [];
cfg.method = 'dfi';
cfg.refindx     = 'all';
cfg.dfi.feature = 'feature';
cfg.dfi.lags    = (1:1:fs/2)./fs;
%dfi1 = ft_connectivityanalysis(cfg, data1);
%dfi2 = ft_connectivityanalysis(cfg, data2);
cfg.dfi.precondition = true;
dfi1b = ft_connectivityanalysis(cfg, data1);
dfi2b = ft_connectivityanalysis(cfg, data2);

cfg = [];
cfg.method = 'mi';
cfg.refindx     = 'all';
cfg.mi.feature = 'feature';
cfg.mi.lags    = (1:1:fs/2)./fs;
%ci1 = ft_connectivityanalysis(cfg, data1);
%ci2 = ft_connectivityanalysis(cfg, data2);
cfg.mi.precondition = true;
ci1b = ft_connectivityanalysis(cfg, data1);
ci2b = ft_connectivityanalysis(cfg, data2);

