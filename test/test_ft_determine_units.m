function test_ft_determine_units

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_determine_units

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
function testElec(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elec = [];
elec.elecpos = randn(32,3);
for i=1:32
  elec.elecpos(i,:) = 100*elec.elecpos(i,:) ./ norm(elec.elecpos(i,:));
end
elec.unit = 'mm';

elec_mm = ft_determine_units(rmfield(ft_convert_units(elec, 'mm'), 'unit'));
elec_cm = ft_determine_units(rmfield(ft_convert_units(elec, 'cm'), 'unit'));
elec_m  = ft_determine_units(rmfield(ft_convert_units(elec, 'm'), 'unit'));

assert(strcmp(elec_mm.unit, 'mm'));
assert(strcmp(elec_cm.unit, 'cm'));
assert(strcmp(elec_m.unit, 'm'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testSphere(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

headmodel = [];
headmodel.r = 120;
headmodel.o = [0 0 40];
headmodel.unit = 'mm';

headmodel_mm = ft_determine_units(rmfield(ft_convert_units(headmodel, 'mm'), 'unit'));
headmodel_cm = ft_determine_units(rmfield(ft_convert_units(headmodel, 'cm'), 'unit'));
headmodel_m  = ft_determine_units(rmfield(ft_convert_units(headmodel, 'm'), 'unit'));

assert(strcmp(headmodel_mm.unit, 'mm'));
assert(strcmp(headmodel_cm.unit, 'cm'));
assert(strcmp(headmodel_m.unit, 'm'));
