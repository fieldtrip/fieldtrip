function test_bug1309

% MEM 2gb
% WALLTIME 00:15:00
% DEPENDENCY

% This tests for the reliability of the ft_convert_units function when dealing
% with different input types (vol,sens,etc.)

elec      = ft_read_sens('standard_1020.elc');
headmodel = ft_read_headmodel('standard_bem.mat');

pos = [20 30 40]; % in mm

cfg = [];
cfg.elec = ft_convert_units(elec, 'mm');
cfg.headmodel = ft_convert_units(headmodel, 'mm');
cfg.dip.pos = pos;
cfg.relnoise = 0.1;
cfg.numtrl = 1;
data = ft_dipolesimulation(cfg);

cfg = [];
timelock = ft_timelockanalysis(cfg, data);

cfg = [];
cfg.model = 'regional';
cfg.latency = [0.195,0.215];
cfg.numdipoles = 1;

% the grid resolution remains the same, regardless of the electrode and headmodel units
% the output source model units should be consistent with the ones specified here (and not with elec or headmodel)
cfg.sourcemodel.resolution = 30;
cfg.sourcemodel.unit = 'mm';

elecunit      = {'mm', 'cm', 'm', 'inch'};
headmodelunit = {'mm', 'cm', 'm', 'feet'};

fprintf('=================================================================================\n');
fprintf('first run to warm-up and precompile\n');
fprintf('=================================================================================\n');
cfg.elec      = ft_convert_units(elec,      elecunit{1});
cfg.headmodel = ft_convert_units(headmodel, headmodelunit{1});
firstrun = ft_dipolefitting(cfg, timelock);

for i=1:numel(elecunit)
  for j=1:numel(headmodelunit)
    stopwatch = tic;
    
    fprintf('=================================================================================\n');
    fprintf('electrodes are expressed in %s, the headmodel is expressed in %s\n', elecunit{i}, headmodelunit{j});
    fprintf('=================================================================================\n');
    cfg.elec      = ft_convert_units(elec,      elecunit{i});
    cfg.headmodel = ft_convert_units(headmodel, headmodelunit{j});
    source{i,j} = ft_dipolefitting(cfg, timelock);
    
    % record the time
    runtime(i,j) = toc(stopwatch);
  end
end

% compute the relative times
runtime = runtime/median(runtime(:));
if any(runtime(:)<0.5) || any(runtime(:)>1.5)
  warning('unexpected difference in execution time for dipole fitting, please check the handling of geometrical units');
end

for i=1:numel(elecunit)
  for j=1:numel(headmodelunit)
    if norm(source{i,j}.dip.pos - pos)>2 % mm
      error('difference in dipole fit results, please check the handling of geometrical units');
    end
  end
end
