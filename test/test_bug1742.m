function test_bug1742

% TEST test_bug1742
% TEST qsubcellfun fexec

% this initially failed due to bug 1742
qsubcellfun(@rand, {1, 2, 3}, 'memreq', 1e9, 'timreq', 600);
% whereas this worked
a = qsubcellfun(@rand, {1, 2, 3}, 'memreq', 1e9, 'timreq', 600);

% perform some more elaborate tests, the svd function behaves differently
% depending on the number of output arguments
            qsubcellfun(@svd, {randn(10, 20), randn(10, 20)}, 'memreq', 1e9, 'timreq', 600);
[s]       = qsubcellfun(@svd, {randn(10, 20), randn(10, 20)}, 'memreq', 1e9, 'timreq', 600);
[u, s]    = qsubcellfun(@svd, {randn(10, 20), randn(10, 20)}, 'memreq', 1e9, 'timreq', 600);
[u, s, v] = qsubcellfun(@svd, {randn(10, 20), randn(10, 20)}, 'memreq', 1e9, 'timreq', 600);
