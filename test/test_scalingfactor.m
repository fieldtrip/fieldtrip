function test_scalingfactor

% MEM 1gb
% WALLTIME 00:03:02

% TEST test_scalingfactor
% TEST scalingfactor ft_convert_units

p = fileparts(which('ft_defaults'));
cd(fullfile(p, 'utilities', 'private'));

assert(scalingfactor('m', 'mm') == 1000);
assert(scalingfactor('mm', 'm') == 0.001);
assert(scalingfactor('T/cm', 'fT/m') == 1e17);

