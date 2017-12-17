function test_bug2888

% WALLTIME 00:10:00
% MEM 1gb

% TEST ft_filetype ft_read_header dataset2files

filename = {
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug2888/0001') % top level directory
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug2888/0001/header')
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug2888/0001/header.txt')
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug2888/0001/samples')
  dccnpath('/home/common/matlab/fieldtrip/data/test/bug2888/0001/events')
  };

for i=1:length(filename)
  cfg = [];
  cfg.dataset = filename{i};
  data = ft_preprocessing(cfg);
  
  hdr = ft_read_header(filename{i});
  dat = ft_read_data(filename{i});
  evt = ft_read_event(filename{i});
end

