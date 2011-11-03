function test_bug1112

% TEST test_bug1112
% TEST ft_convert_units ft_estimate_units

load test_bug1112
sens = ft_convert_units(sens);

if ~strcmp(sens.unit, 'cm')
  error('the estimated units are incorrect');
end
