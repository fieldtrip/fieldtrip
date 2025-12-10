function test_issue2075

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY nanmean
% DATA no

% A user reported an issue related to nanmean, using MATLAB 2022a.
% Operating system unknown. This script tries to reproduce it. On our end,
% running on Linux we cannot reproduce an error (the exact error is as of
% yet still underspecified). The purpose of this function is to document
% it, being aware that the test functions are not necessarily executed in
% their automatic overnight batch in the reported matlab version

[ftver, ftpath] = ft_version;

global ft_default
ft_default.toolbox.signal = 'compat';
ft_default.toolbox.stats= 'compat';
restoredefaultpath;
addpath(ftpath);
ft_defaults;
ft_preproc_bandpassfilter(randn(1,1000), 1000, [230 0.5], [], 'fir');
which nanmean

ft_default.toolbox.signal = 'matlab'; % rather than compat
ft_default.toolbox.stats= 'matlab';
restoredefaultpath;
addpath(ftpath);
ft_defaults;
ft_preproc_bandpassfilter(randn(1,1000), 1000, [230 0.5], [], 'fir');
which nanmean
