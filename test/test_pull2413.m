function test_pull2413

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_event
% DATA private

%%

datadir = dccnpath('/project/3031000.02/test/original/eeg/gdf/');

d = dir(fullfile(datadir, '*.gdf'));
for k = 1:numel(d)
  event = ft_read_event(fullfile(d(k).folder, d(k).name));
end

