function test_ft_convert_units

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_convert_units

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
elec.chanpos = [10 10 10];
elec.elecpos = [10 10 10];
elec.unit = 'mm';

elec_um = ft_convert_units(elec, 'um');
elec_mm = ft_convert_units(elec, 'mm');
elec_cm = ft_convert_units(elec, 'cm');
elec_dm = ft_convert_units(elec, 'dm');
elec_m  = ft_convert_units(elec, 'm');
elec_hm = ft_convert_units(elec, 'hm');
elec_km = ft_convert_units(elec, 'km');

assert(all(round(elec_um.chanpos ./ elec.chanpos, 0) == 1e+3)); % micrometer
assert(all(round(elec_mm.chanpos ./ elec.chanpos, 0) == 1e-0));
assert(all(round(elec_cm.chanpos ./ elec.chanpos, 1) == 1e-1));
assert(all(round(elec_dm.chanpos ./ elec.chanpos, 2) == 1e-2));
assert(all(round(elec_m.chanpos  ./ elec.chanpos, 3) == 1e-3));
assert(all(round(elec_hm.chanpos ./ elec.chanpos, 5) == 1e-5)); % hectometer
assert(all(round(elec_km.chanpos ./ elec.chanpos, 6) == 1e-6));

assert(all(round(elec_um.elecpos ./ elec.elecpos, 0) == 1e+3));
assert(all(round(elec_mm.elecpos ./ elec.elecpos, 0) == 1e-0));
assert(all(round(elec_cm.elecpos ./ elec.elecpos, 1) == 1e-1));
assert(all(round(elec_dm.elecpos ./ elec.elecpos, 2) == 1e-2));
assert(all(round(elec_m.elecpos  ./ elec.elecpos, 3) == 1e-3));
assert(all(round(elec_hm.elecpos ./ elec.elecpos, 5) == 1e-5));
assert(all(round(elec_km.elecpos ./ elec.elecpos, 6) == 1e-6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function testSphere(testCase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

headmodel = [];
headmodel.r = 120;
headmodel.o = [0 0 40];
headmodel.unit = 'mm';

headmodel_um = ft_convert_units(headmodel, 'um');
headmodel_mm = ft_convert_units(headmodel, 'mm');
headmodel_cm = ft_convert_units(headmodel, 'cm');
headmodel_dm = ft_convert_units(headmodel, 'dm');
headmodel_m  = ft_convert_units(headmodel, 'm');
headmodel_hm = ft_convert_units(headmodel, 'hm');
headmodel_km = ft_convert_units(headmodel, 'km');

assert(all(round(headmodel_um.r ./ headmodel.r, 0) == 1e+3)); % micrometer
assert(all(round(headmodel_mm.r ./ headmodel.r, 0) == 1e-0));
assert(all(round(headmodel_cm.r ./ headmodel.r, 1) == 1e-1));
assert(all(round(headmodel_dm.r ./ headmodel.r, 2) == 1e-2));
assert(all(round(headmodel_m.r  ./ headmodel.r, 3) == 1e-3));
assert(all(round(headmodel_hm.r ./ headmodel.r, 5) == 1e-5)); % hectometer
assert(all(round(headmodel_km.r ./ headmodel.r, 6) == 1e-6));

assert(all(round(headmodel_um.o(3) ./ headmodel.o(3), 0) == 1e+3));
assert(all(round(headmodel_mm.o(3) ./ headmodel.o(3), 0) == 1e-0));
assert(all(round(headmodel_cm.o(3) ./ headmodel.o(3), 1) == 1e-1));
assert(all(round(headmodel_dm.o(3) ./ headmodel.o(3), 2) == 1e-2));
assert(all(round(headmodel_m.o(3)  ./ headmodel.o(3), 3) == 1e-3));
assert(all(round(headmodel_hm.o(3) ./ headmodel.o(3), 5) == 1e-5));
assert(all(round(headmodel_km.o(3) ./ headmodel.o(3), 6) == 1e-6));
