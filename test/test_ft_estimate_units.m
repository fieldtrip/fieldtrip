function test_ft_estimate_units

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_estimate_units

if nargout
  % assume that this is called by RUNTESTS
  tests = functiontests(localfunctions);
else
  % assume that this is called from the command line
  func = localfunctions;
  for i=1:numel(func)
    fprintf('evaluating %s\n', func2str(func{i}));
    feval(func{i});
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testUnits(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(strcmp(ft_estimate_units(80),  'mm'));
assert(strcmp(ft_estimate_units(8),   'cm'));
assert(strcmp(ft_estimate_units(0.08), 'm'));

assert(strcmp(ft_estimate_units(100),  'mm'));
assert(strcmp(ft_estimate_units(10),   'cm'));
assert(strcmp(ft_estimate_units(0.10), 'm'));

assert(strcmp(ft_estimate_units(120),  'mm'));
assert(strcmp(ft_estimate_units(12),   'cm'));
assert(strcmp(ft_estimate_units(0.12), 'm'));
