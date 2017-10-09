function test_bug3361

% MEM 2gb
% WALLTIME 00:10:00

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3361.mat'));

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
cfg.highlightchannel = 'all';

close all
ft_topoplotTFR(cfg, ga)
