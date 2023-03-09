function test_pull2183

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_headshape ft_read_headmodel

simnibsdir = dccnpath('/home/common/matlab/fieldtrip/data/test/original/simnibs/');
d = dir(fullfile(simnibsdir, '*.msh'));
for k = 1:numel(d)
  if strcmp(d(k).name, 'sphere3_v4_1_ascii.msh')
    continue;
  elseif strcmp(d(k).name, 'sphere3_v4_1_binary.msh')
    continue;
  end
  test = ft_read_headshape(fullfile(d(k).folder, d(k).name));
end

pwdir = pwd;
cd(simnibsdir);
cd('v4');
hm = ft_read_headmodel('ernie.msh', 'fileformat', 'simnibs');
hm = ft_read_headmodel('ernie.msh', 'fileformat', 'simnibs', 'meshtype', 'surface');
cd(simnibsdir);
cd('v3');
hm = ft_read_headmodel('ernie_v3.msh', 'fileformat', 'simnibs');
hm = ft_read_headmodel('ernie_v3.msh', 'fileformat', 'simnibs', 'meshtype', 'surface');

cd(pwdir);