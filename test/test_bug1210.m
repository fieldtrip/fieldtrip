function test_bug1210

% this is a test script that relates to http://bugzilla.fcdonders.nl/show_bug.cgi?id=1210

% TEST test_bug1210
% TEST ft_datatype_sens

load test_bug1210

% the input is old-style
assert(isfield(D.grad, 'pnt'));
assert(isfield(D.grad, 'ori'));

gradnew = ft_datatype_sens(D.grad);

% the output should be new-style
assert(isfield(gradnew, 'chanpos'));
assert(isfield(gradnew, 'chanori'));
assert(isfield(gradnew, 'coilpos'));
assert(isfield(gradnew, 'coilori'));
