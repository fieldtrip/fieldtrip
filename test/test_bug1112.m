function test_bug1112

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_convert_units ft_estimate_units

cd(dccnpath('/home/common/matlab/fieldtrip/data/test'))
load bug1112

sens = ft_convert_units(sens);

if ~strcmp(sens.unit, 'cm')
  error('the estimated units are incorrect');
end
