function test_bug2231

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_header ft_preprocessing read_bti_m4d
% DATA private

% Bug reported by Christian Wienbruch, about the functionality of reading 4D-data
% the old-fashioned way being broken

datadir = dccnpath('/home/common/matlab/fieldtrip/data/test/bug2231/');
d = dir(fullfile(datadir, '*,s'));
filename = fullfile(datadir, d(1).name);

cfg = [];
cfg.dataset = filename;
data = ft_preprocessing(cfg);
hdr  = ft_read_header(cfg.dataset);

% if we get here, at least there has been no crash
fprintf('Reading 4D data the old fashioned way did not lead to a crash\n');
