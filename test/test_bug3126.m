function test_bug3126

% WALLTIME 00:20:00
% MEM 4gb

filename = dccnpath('/home/common/matlab/fieldtrip/data/test/original/eeg/neurone/2016-02-11T152726');

%% low level functions

hdr = ft_read_header(filename);
dat = ft_read_data(filename);
evt = ft_read_event(filename);

%% read the continuous data as one segment

cfg = [];
cfg.dataset = filename;
data1 = ft_preprocessing(cfg);

%% read the first 100 seconds as trials

cfg = [];
cfg.dataset = filename;
cfg.trl(:,1) = (0:99)*1000 + 1;
cfg.trl(:,2) = (0:99)*1000 + 1000;
cfg.trl(:,3) = 0;
data2 = ft_preprocessing(cfg);


%% browse the data

if false
  % this section should not run automatically
  cfg = [];
  cfg.viewmode = 'vertical';
  cfg.dataset = filename;
  ft_databrowser(cfg);
end
