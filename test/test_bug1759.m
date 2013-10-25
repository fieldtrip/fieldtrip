function test_bug1759

% MEM 1gb
% WALLTIME 00:03:01

% TEST test_bug1759

% Sparse matrix multplication results in slightly different results than nonsparse
% http://bugzilla.fcdonders.nl/show_bug.cgi?id=1759

cd(dccnfilename('/home/common/matlab/fieldtrip/data/test'))
load bug1759.mat

tmp{1} = tra*dat;
tmp{2} = sparse(tra)*dat;

assert(identical(tmp{1}, tmp{2}, 'reltol', 0.0001));

