function test_bug1910

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_datatype_sens

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug1910.mat'));

% update the old gradiometer structure
gradnew = ft_datatype_sens(grad);

disp(gradnew);

assert(all(~isnan(gradnew.chanpos(:))))
