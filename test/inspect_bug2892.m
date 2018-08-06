function inspect_bug2892

%% [p, f, x] = fileparts(mfilename('fullpath'));
p = tempname; 
mkdir(p);
cd(p)
fprintf('temporary results are stored in %s\n', p);

%% do not use local variables, only data on disk

cfg = [];
cfg.numtrl = 100;
cfg.outputfile = 'data1.mat';
ft_freqsimulation(cfg);

cfg.outputfile = 'data2.mat';
ft_freqsimulation(cfg); % different random noise

cfg = [];
cfg.method = 'wavelet';
cfg.toi = 0:0.01:1;
cfg.channel = 1;
cfg.keeptrials = 'yes';
cfg.inputfile  = 'data1.mat';
cfg.outputfile = 'freq1.mat';
ft_freqanalysis(cfg);

cfg.inputfile  = 'data2.mat';
cfg.outputfile = 'freq2.mat';
ft_freqanalysis(cfg);


cfg = [];
cfg.operation = 'multiply';
cfg.parameter = 'powspctrm';
cfg.scalar = 0.9;
cfg.inputfile  = 'freq1.mat';
cfg.outputfile = 'cond1.mat';
ft_math(cfg);

cfg.scalar = 1.1;
cfg.inputfile  = 'freq2.mat';
cfg.outputfile = 'cond2.mat';
ft_math(cfg);

cfg = [];
cfg.statistic = 'indepsamplesT';
cfg.method = 'analytic';
cfg.ivar = 1;
cfg.design = [1*ones(1,100) 2*ones(1,100)];
cfg.inputfile  = {'freq1.mat', 'freq2.mat'};
cfg.outputfile = 'stat.mat';
ft_freqstatistics(cfg);

cfg = [];
cfg.inputfile = 'stat.mat';
ft_analysispipeline(cfg);

