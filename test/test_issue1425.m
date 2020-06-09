function test_issue1425(testfilename)

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_preprocessing ft_read_data

if nargin==0
  d = dir(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1425/*.edf'));
  testfilename = fullfile(d.folder, d.name);
end
hdr   = ft_read_header(fullfile(testfilename));
event = ft_read_event(fullfile(testfilename));

