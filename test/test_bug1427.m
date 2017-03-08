function test_bug1427

% MEM 2gb
% WALLTIME 00:10:00

% TEST dataset2files ft_read_header ft_read_data

% the Long64ChannelWithEvents fails because it consists of muliple segments
% filepath = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1427/Long64ChannelWithEvents.mff');

filepath = dccnpath('/home/common/matlab/fieldtrip/data/test/bug1427/NS500Sine6Hz.mff');

filename = dir(filepath);
filename = {filename(~[filename.isdir]).name};

% determine the type for any of the files
typ = cell(size(filename));
for i=1:length(filename)
  typ{i} = ft_filetype(fullfile(filepath, filename{i}));
end

sel = find(strcmp(typ, 'unknown'));
if ~isempty(sel)
  error('the filetype detection failed on %s',  fullfile(filepath, filename{sel(1)}))
end

% read the header for any of the files
hdr = cell(size(filename));
for i=1:length(filename)
  disp(filename{i})
  hdr{i} = ft_read_header(fullfile(filepath, filename{i}));
end

% read some data for any of the files
dat = cell(size(filename));
for i=1:length(filename)
  disp(filename{i})
  dat{i} = ft_read_data(fullfile(filepath, filename{i}), 'begsample', 1, 'endsample', 1000);
end

