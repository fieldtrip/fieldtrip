function test_bug1210

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_datatype_sens
% DATA private

cd(dccnpath('/project/3031000.02/test'))
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
