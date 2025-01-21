function test_bug1112

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_convert_units ft_estimate_units
% DATA private

cd(dccnpath('/project/3031000.02/test'))
load bug1112

sens = ft_determine_units(sens);

if ~strcmp(sens.unit, 'cm')
  error('the estimated units are incorrect');
end
