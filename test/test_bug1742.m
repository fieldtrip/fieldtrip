function test_bug1742

% MEM 1500mb
% WALLTIME 00:10:00

% TEST qsubcellfun fexec

[ftver, ftpath] = ft_version;
addpath(fullfile(ftpath,'qsub'));

% this initially failed due to bug 1742
qsubcellfun(@rand, {1, 2, 3}, 'memreq', 1e9, 'timreq', 600, 'backend', 'system');
% whereas this worked
a = qsubcellfun(@rand, {1, 2, 3}, 'memreq', 1e9, 'timreq', 600, 'backend', 'system');

% perform some more elaborate tests, the svd function behaves differently
% depending on the number of output arguments
            qsubcellfun(@svd, {randn(10, 20), randn(10, 20)}, 'memreq', 1e9, 'timreq', 600, 'backend', 'local');
[s]       = qsubcellfun(@svd, {randn(10, 20), randn(10, 20)}, 'memreq', 1e9, 'timreq', 600, 'backend', 'local');
[u, s]    = qsubcellfun(@svd, {randn(10, 20), randn(10, 20)}, 'memreq', 1e9, 'timreq', 600, 'backend', 'local');
[u, s, v] = qsubcellfun(@svd, {randn(10, 20), randn(10, 20)}, 'memreq', 1e9, 'timreq', 600, 'backend', 'local');
