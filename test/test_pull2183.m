function test_pull2183

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_headshape ft_read_headmodel
% PRIVATEDATA

%%

mshdir = dccnpath('/home/common/matlab/fieldtrip/data/test/pull2183/');

d = dir(fullfile(mshdir, '*.msh'));
for k = 1:numel(d)
  hs = ft_read_headshape(fullfile(d(k).folder, d(k).name));
end

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/simnibs/v4/ernie.msh');
hm = ft_read_headmodel(filename, 'fileformat', 'simnibs');
hm = ft_read_headmodel(filename, 'fileformat', 'simnibs', 'meshtype', 'surface');

%%

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/simnibs/v3/ernie_v3.msh');
hm = ft_read_headmodel(filename, 'fileformat', 'simnibs');
hm = ft_read_headmodel(filename, 'fileformat', 'simnibs', 'meshtype', 'surface');

