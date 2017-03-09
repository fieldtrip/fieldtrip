function test_bug1210

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_datatype_sens

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1210.mat

% the input is old-style
assert(isfield(D.grad, 'pnt'));
assert(isfield(D.grad, 'ori'));

gradnew = ft_datatype_sens(D.grad);

% the output should be new-style
assert(isfield(gradnew, 'chanpos'));
assert(isfield(gradnew, 'chanori'));
assert(isfield(gradnew, 'coilpos'));
assert(isfield(gradnew, 'coilori'));
