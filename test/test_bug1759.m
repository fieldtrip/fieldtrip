function test_bug1759

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY


% Sparse matrix multplication results in slightly different results than nonsparse
% http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=1759

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1759.mat

tmp{1} = tra*dat;
tmp{2} = sparse(tra)*dat;

assert(isalmostequal(tmp{1}, tmp{2}, 'reltol', 0.0001));

