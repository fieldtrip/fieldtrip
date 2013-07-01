function test_bug1243

% TEST test_bug1243
% TEST ft_topoplotIC

cd(dccnfilename('/home/common/matlab/fieldtrip/data/test'))
load bug1243.mat

figure
for i=1:9
  subplot(3,3,i);
  cfg = [];
  cfg.component = i;
  cfg.layout = 'EEG1020.lay';
  ft_topoplotIC(cfg, comp);
end
