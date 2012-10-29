function test_scalingfactor

% TEST test_scalingfactor
% TEST scalingfactor ft_convert_units

assert(scalingfactor('m', 'mm') == 1000);
assert(scalingfactor('mm', 'm') == 0.001);
assert(scalingfactor('T/cm', 'fT/m') == 1e17);

