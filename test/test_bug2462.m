function test_bug2462

% MEM 1500mb
% WALLTIME 00:10:00
% test_bug2462

homedir = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2462/');
datasets = {
  'scan1_Filters_125HzLP-ascii-multiplexed.dat'
  'scan1_Filters_125HzLP-ascii-vectorized.dat'
  'scan1_Filters_125HzLP-binary-multiplexed.dat'
  'scan1_Filters_125HzLP-binary-vectorized.dat'
};

data = cell(numel(datasets),1);

for i=1:numel(datasets)
  cfg = [];
  cfg.dataset = fullfile(homedir, datasets{i});
  
  data{i} = ft_preprocessing(cfg);
end;


tolerance = 1;
% Fixme - create and load different brainvision datasets.
for i=1:numel(data)
  for j=1:numel(data)
    assert(sum(abs(data{i}.trial{1}(:)-data{j}.trial{1}(:)))<tolerance);
  end
end;

end


