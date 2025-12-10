function test_pull1826

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_read_headshape
% DATA private

%%

datadir = dccnpath('/project/3031000.02/test/original/headshape/structuresensor/');
folders = strsplit(ls(datadir));

for k = 1:(numel(folders)-1)
  % try if it loads without error
  tmp = ft_read_headshape(fullfile(datadir, folders{k}, 'Model.obj'));
end
