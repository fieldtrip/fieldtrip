function inspect_qsubcellfun5

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY

% this test script pertains to
% http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1255

%%

y = qsubcellfun(@inspect_qsubcellfun5_fn, {1, 2, 3}, 'memreq', 0, 'timreq', 500)
y = qsubcellfun(@inspect_qsubcellfun5_fn, {1, 2, 3}, 'memreq', 0, 'timreq', 500, 'stack', 3)

%%

fc = qsubcompile(@inspect_qsubcellfun5_fn)
y = qsubcellfun(fc, {1, 2, 3}, 'memreq', 0, 'timreq', 500)

%%

y = qsubcellfun(@inspect_qsubcellfun5_fn, {1, 2, 3}, 'memreq', 0, 'timreq', 500, 'compile', true, 'stack', 1)
