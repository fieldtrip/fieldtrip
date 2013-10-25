function test_bug1910

% MEM 1gb
% WALLTIME 00:03:01

% TEST test_bug1910
% TEST ft_datatype_sens

load(dccnfilename('/home/common/matlab/fieldtrip/data/test/bug1910.mat'));

% update the old gradiometer structure
gradnew = ft_datatype_sens(grad);

disp(gradnew);

assert(all(~isnan(gradnew.chanpos(:))))
