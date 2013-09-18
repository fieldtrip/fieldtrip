function test_bug1984_2187

% TEST test_bug1984_2187
% TEST ft_appendfreq ft_freqgrandaverage ft_freqstatistics ft_prepare_neighbours

cd(dccnfilename('/home/common/matlab/fieldtrip/data/test'));
load bug1984_2187.mat

%% first: test rpt concatenation
cfg=[];
cfg.parameter='powspctrm';
freq = ft_appendfreq(cfg,s,s,s);
if ~isfield(freq,'freq');
  error('freq field is not appeneded: see bugs 1984 and 2187');
end

%% first: check functionality when using ft_freqstatistis using ft_freqgrandaverage before
cfg = [];
cfg.keepindividual = 'yes';
c1 = ft_freqgrandaverage(cfg,cond1{:});
c2 = ft_freqgrandaverage(cfg,cond2{:});

cfg                  = [];
cfg.channel          = {'MEG'};
cfg.latency          = [1.8 2.1];
cfg.frequency        = [60 90];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum'; 
cfg.correcttail      = 'alpha';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.clusterthreshold = 'nonparametric_common';
cfg.minnbchan        = 2;
cfg.numrandomization = 100;
cfg_neighb.method    = 'template';
cfg_neighb.template  = 'CTF275_neighb.mat';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb);

subj = 5;
design = zeros(2,2*subj);
for i = 1:subj
  design(1,i) = i;
end
for i = 1:subj
  design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design           = design;
cfg.ivar             = 2;
cfg.uvar             = 1;
cfg.avgoverfreq      = 'no';
cfg.avgovertime      = 'yes';
cfg.parameter        = 'powspctrm';

stat = ft_freqstatistics(cfg, c1, c2);
if ~isfield(stat,'freq');
  error('freq field is not appeneded: might be related with ft_appendfreq bugs 1984 and 2187');
end

%% second: check functionality when using ft_freqstatistis using the varargin
stat = ft_freqstatistics(cfg, cond1{:}, cond2{:});
if ~isfield(stat,'freq');
  error('freq field is not appeneded: might be related with ft_appendfreq bugs 1984 and 2187');
end

