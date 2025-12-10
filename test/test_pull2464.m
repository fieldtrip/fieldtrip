function test_pull2464

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_read_headshape
% DATA private

%%

datadir = dccnpath('/project/3031000.02/test/');
fname   = fullfile(datadir, 'pull2464.fif');
hs      = ft_read_headshape(fname, 'format', 'mne_source');

