function test_scalingfactor

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_scalingfactor ft_convert_units

% since the function to test is in a private directory, we explicitely have to cd into that directory
[ftver, ftpath] = ft_version;
cd(fullfile(ftpath, 'utilities', 'private'));

assert(ft_scalingfactor('m', 'mm') == 1000);
assert(ft_scalingfactor('mm', 'm') == 0.001);
assert(ft_scalingfactor('T/cm', 'fT/m') == 1e17);

