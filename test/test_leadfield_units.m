% function test_leadfield_units

% TEST test_leadfield_units
% TEST ft_compute_leadfield FT_PREPARE_VOL_SENS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sens          = [];
sens.coilpos  = [0.01 0 0.01; -0.01 0 0.01];
sens.coilori  = [0 0 1; 0 0 1];
sens.tra      = [1 -1] / 0.02; % divide by the baseline
sens.label    = {'M1'};
sens.chantype = {'megplanar'};
sens.chanunit = {'T/m'};
sens.unit     = 'm';
sens.type     = 'meg';

if false
  sens          = [];
  sens.coilpos  = [0 0 0.01];
  sens.coilori  = [0 0 1];
  sens.chanpos  = [0 0 0.01];
  sens.chanori  = [0 0 1];
  sens.tra      = 1;
  sens.label    = {'M1'};
  sens.chantype = {'megmag'};
  sens.chanunit = {'T'};
  sens.unit     = 'm';
  sens.type     = 'meg';
end

vol      = [];
vol.type = 'infinite';
vol.unit = 'm';

pos = [0 0 -0.01];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos_m  = pos;
pos_cm = pos*100;
pos_mm = pos*1000;

vol_m  = ft_convert_units(vol, 'm');
vol_cm = ft_convert_units(vol, 'cm');
vol_mm = ft_convert_units(vol, 'mm');

sens_m  = ft_convert_units(sens, 'm');
sens_cm = ft_convert_units(sens, 'cm');
sens_mm = ft_convert_units(sens, 'mm');

if strcmp(sens.chantype{1}, 'megplanar')
  % this should be scaled with the geometrical units
  assert(~isequal(sens_m.tra, sens_cm.tra)); 
  assert(~isequal(sens_m.tra, sens_cm.tra)); 
end

arbitrary.lf_m  = ft_compute_leadfield(pos_m , sens_m , vol_m , 'unit', 'arbitrary');
arbitrary.lf_cm = ft_compute_leadfield(pos_cm, sens_cm, vol_cm, 'unit', 'arbitrary');
arbitrary.lf_mm = ft_compute_leadfield(pos_mm, sens_mm, vol_mm, 'unit', 'arbitrary');

si.lf_m  = ft_compute_leadfield(pos_m , sens_m , vol_m , 'unit', 'si');
si.lf_cm = ft_compute_leadfield(pos_cm, sens_cm, vol_cm, 'unit', 'si');
si.lf_mm = ft_compute_leadfield(pos_mm, sens_mm, vol_mm, 'unit', 'si');

