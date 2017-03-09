function test_bug1315

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_checkdata ft_prepare_neighbours ft_megplanar ft_combineplanar

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1315.mat'))

% neighbours
cfg = [];
cfg.method = 'template';
cfg.layout = 'CTF275.lay';
neighbours = ft_prepare_neighbours(cfg, data);
%% producing the bug
cfg = [];
cfg.neighbours = neighbours;
erf_planar = ft_megplanar(cfg, data);

erf_combined = ft_combineplanar([], erf_planar);

if size(erf_combined.time, 2) ~= size(erf_combined.trial, 3) ...
   || any(erf_combined.time ~= data.time)
  error('Time axis screwed up (again)');
end

  
