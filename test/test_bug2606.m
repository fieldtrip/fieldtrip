function test_bug2606

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_timelockanalysis

%% generate some data

% get CTF275 labels
lay = ft_prepare_layout(struct('layout', 'CTF275.lay'));
labels = lay.label(1:275);

nsubj = 20;
condA = cell(nsubj,1);
condB = cell(nsubj,1);

for k = 1:nsubj
  tl = [];
  tl.label = labels;
  tl.time = 0.001:0.001:1;
  tl.avg = randn(numel(labels), numel(tl.time));
  tl.dimord = 'chan_time';
  
  condA{k} = tl;
  
  tl.avg = randn(numel(labels), numel(tl.time)) + 1;
  condB{k} = tl;
end

%% select these sensors to average over

sensors_lefttemp = [
   104
    99
   100
   108
   101
   106
   109
   102
   110
   107
   ];

%% call timelockstatistics

cfg = [];
cfg.avgovertime = 'no';
cfg.parameter = 'avg';
cfg.method = 'montecarlo'; % see ft_statistics_montecarlo

cfg.correctm = 'cluster'; % multiple-comparison correction
cfg.statistic = 'depsamplesT';
cfg.tail = 0; % two-sided test
cfg.alpha = 0.05;
cfg.correcttail = 'alpha'; % bonferroni correction

cfg.clusteralpha = 0.05; 
cfg.clusterstatistic = 'maxsum';
% load('/home/predatt/heisol/3018012.05_heisol/meg/gavg/neighbours.mat');

c2 = [];
c2.method = 'template';
c2.template = 'ctf275_neighb.mat';
neighbours = ft_prepare_neighbours(c2);

cfg.neighbours = neighbours;

design = zeros(2,2*nsubj);
for i = 1:nsubj
design(1,i) = i;
end
for i = 1:nsubj
design(1,nsubj+i) = i;
end
design(2,1:nsubj)        = 1;
design(2,nsubj+1:2*nsubj) = 2;

cfg.design = design;
cfg.uvar  = 1; % "unit variable" (subject)
cfg.ivar  = 2; % independent variable (condition)

cfg.numrandomization = 1000;

cfg.latency = [0.1 0.4];
cfg.channel = sensors_lefttemp;
cfg.avgoverchan = 'yes';

Test = ft_timelockstatistics(cfg, condA{:}, condB{:}); 
