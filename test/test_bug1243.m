function test_bug1243

load test_bug1243;
figure
for i=1:9
  subplot(3,3,i);
  cfg = [];
  cfg.component = i;
  cfg.layout = 'EEG1020.lay';
  ft_topoplotIC(cfg, comp);
end
