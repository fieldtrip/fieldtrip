function test_bug2231

% MEM 1500mb
% WALLTIME 00:10:00

% TEST: ft_read_header
% TEST: ft_preprocessing
% TEST: read_bti_m4d

% Bug reported by Christian Wienbruch, about the functionality of reading 4D-data
% the old-fashioned way being broken

d = dir(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2231/*,s'));
filename = fullfile(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2231/'), d(1).name);

cfg = [];
cfg.dataset = filename;
data = ft_preprocessing(cfg);
hdr  = ft_read_header(cfg.dataset);

% if we get here, at least there has been no crash
fprintf('Reading 4D data the old fashioned way did not lead to a crash\n');
