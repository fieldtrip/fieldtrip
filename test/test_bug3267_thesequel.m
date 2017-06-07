function test_bug3267_thesequel

% MEM 1500mb
% WALLTIME 00:10:00


%% situation 1: without sampleinfos (works fine)
data1a.label = {'1';'2';'3'};
for trl = 1:5
  data1a.trial{1,trl} = randn(3,3);
  data1a.time{1,trl} = [1 2 3];
end
data1a.trialinfo(:,1) = 1:5;

data1b = data1a;
data1b.trialinfo(:,1) = 6:10;

% append
data1 = ft_appenddata([], data1a, data1b);

% this works fine
cfg             = [];
cfg.length      = 1;
cfg.overlap     = 0;
redef1 = ft_redefinetrial(cfg, data1);

% this should be 1 1 1 2 2 2 ... 10 10 10
assert(all(abs(diff(redef1.trialinfo))<2));


%% situation 2: with identical sampleinfos (produces incorrect trialinfo)
data2a = data1a;
data2a.sampleinfo = [1 3; 4 6; 8 10; 11 13; 16 18];
data2b = data1b;
data2b.sampleinfo = data2a.sampleinfo;

% append
data2 = ft_appenddata([], data2a, data2b);

% tacitly fails because of intertwined sampleinfos (at line 263)
cfg             = [];
cfg.length      = 1;
cfg.overlap     = 0;
redef2 = ft_redefinetrial(cfg, data2);

% this should be 1 1 1 2 2 2 ... 10 10 10
assert(all(abs(diff(redef2.trialinfo))<2));


%% situation 3: with intertwined and partly overlapping sampleinfos (produces incorrect trialinfo, as in situation 2)
data3a = data1a;
data3a.sampleinfo = [1 3; 4 6; 8 10; 11 13; 16 18];
data3b = data1b;
data3b.sampleinfo = [2 4; 5 7; 9 11; 12 14; 15 17];

% append
data3 = ft_appenddata([], data3a, data3b);

% tacitly fails because of intertwined sampleinfos (at line 263)
cfg             = [];
cfg.length      = 1;
cfg.overlap     = 0;
redef3 = ft_redefinetrial(cfg, data3);

% this should be 1 1 1 2 2 2 ... 10 10 10
assert(all(abs(diff(redef3.trialinfo))<2));

