function test_pull2582

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_data ft_read_event ft_read_header 
% DATA no

datadir = dccnpath('/project/3031000.02/test/original/eeg/brainvision');

ft_hastoolbox('bvrf_reader',1);

cd(fullfile(datadir, 'Noise'));
dataset = '001.bvrd';

cfg = [];
cfg.dataset = dataset;
data = ft_preprocessing(cfg);
assert(numel(data.time{1})==112930); % this catches the buggy original code from brainvision

hdr   = ft_read_header('001.bvrh');
event = ft_read_event('001.bvrm');

cd(fullfile(datadir, 'DemoR2_VisualOddball'));
dataset = 'DemoR2_VisualOddball_actiCHamp_slim_42_01.bvrd';
cfg.dataset = dataset;
data = ft_preprocessing(cfg);

dataset = 'DemoR2_VisualOddball_actiCHamp_slim_42_02.bvrd';
cfg.dataset = dataset;
data = ft_preprocessing(cfg);

hdr = ft_read_header('DemoR2_VisualOddball_actiCHamp_slim_42_01.bvrh');
event = ft_read_event('DemoR2_VisualOddball_actiCHamp_slim_42_01.bvrm');
