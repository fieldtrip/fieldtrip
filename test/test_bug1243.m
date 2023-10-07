function test_bug1243

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_topoplotIC
% DATA private

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1243.mat'));

figure
for i=1:9
  cfg = [];
  cfg.component = i;
  cfg.figure = subplot(3,3,i);
  cfg.layout = 'EEG1020.lay';
  ft_topoplotIC(cfg, comp);
end
