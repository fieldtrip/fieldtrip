function test_bug3361

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY

% I noticed an error in MATLAB 2017b which would try to evaluate "ga" as a
% function, even though it was read from the file. The solution is to read
% it explicitly.
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3361.mat'), 'cfg');
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3361.mat'), 'ga');

%%

close all
ft_topoplotTFR(cfg, ga)

%%

cfg.highlight = {'on', 'on'};
cfg.highlightchannel = {
  {'MLO11'
  'MLO12'
  'MLO13'
  'MLO14'
  'MLO21'
  'MLO22'
  'MLO23'
  'MLO24'
  'MLO31'
  'MLO32'
  'MLO34'
  'MLO41'
  'MLO42'
  'MLO43'
  'MLO44'
  'MLO51'
  'MLO52'
  'MLO53'}
  {'MRO11'
  'MRO12'
  'MRO13'
  'MRO14'
  'MRO21'
  'MRO22'
  'MRO23'
  'MRO24'
  'MRO31'
  'MRO32'
  'MRO33'
  'MRO34'
  'MRO41'
  'MRO42'
  'MRO43'
  'MRO44'
  'MRO51'
  'MRO52'
  'MRO53'}};
cfg.highlightsymbol = {'o', '*'};

close all
ft_topoplotTFR(cfg, ga)

%%

cfg.highlight = 'on';
cfg.highlightsymbol = '+';
cfg.highlightsize = 30;
cfg.highlightchannel = 'all';

close all
ft_topoplotTFR(cfg, ga)


%%
% I detected another problem with clusterplot and highlighting

clear all
close all

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3361b.mat'), 'cfg');
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3361b.mat'), 'stat');

% fake the first 5 positive clusters such that they show up
stat.posclusters(1).prob = 0.01 - 0.0001;
stat.posclusters(2).prob = 0.05 - 0.0001;
stat.posclusters(3).prob = 0.10 - 0.0001;
stat.posclusters(4).prob = 0.20 - 0.0001;
stat.posclusters(5).prob = 0.30 - 0.0001;
% ensure that the negative clusters don't show up
stat.negclusters(1).prob = 1;
stat.negclusters(2).prob = 1;
stat.negclusters(3).prob = 1;
stat.negclusters(4).prob = 1;
stat.negclusters(5).prob = 1;

cfg.alpha = 0.3;
ft_clusterplot(cfg, stat);
