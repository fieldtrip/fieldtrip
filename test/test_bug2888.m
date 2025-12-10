function test_bug2888

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_filetype ft_read_header dataset2files
% DATA private

filename = {
  dccnpath('/project/3031000.02/test/bug2888/0001') % top level directory
  dccnpath('/project/3031000.02/test/bug2888/0001/header')
  dccnpath('/project/3031000.02/test/bug2888/0001/header.txt')
  dccnpath('/project/3031000.02/test/bug2888/0001/samples')
  dccnpath('/project/3031000.02/test/bug2888/0001/events')
  };

for i=1:length(filename)
  cfg = [];
  cfg.dataset = filename{i};
  data = ft_preprocessing(cfg);
  
  hdr = ft_read_header(filename{i});
  dat = ft_read_data(filename{i});
  evt = ft_read_event(filename{i});
end

