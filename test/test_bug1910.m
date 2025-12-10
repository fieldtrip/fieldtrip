function test_bug1910

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_datatype_sens
% DATA private

load(dccnpath('/project/3031000.02/test/bug1910.mat'));

% update the old gradiometer structure
gradnew = ft_datatype_sens(grad);

disp(gradnew);

assert(all(~isnan(gradnew.chanpos(:))))
