function test_bug1309

% MEM 1500mb
% WALLTIME 00:10:00


% This tests for the reliability of the ft_convert_units function when dealing
% with different input types (vol,sens,etc.)

global ft_default;
ft_default.feedback = 'no';

% artificially create timelock data
timelock_data.fsamepl = 500;
timelock_data.dimord = 'chan_time';
timelock_data.time = zeros(1, 500);
timelock_data.avg = randn(31, 500);
for i = 1:31
  timelock_data.label{i} = sprintf('label%d', i);
end

elec      = ft_read_sens('standard_1020.elc');
headmodel = ft_read_vol('standard_bem.mat');

timelock_data.label = elec.label(1:31);

cfg = [];
cfg.model = 'regional';
cfg.latency = [0.14,0.17];
cfg.numdipoles = 1;

% the grid resolution remains the same, regardless of the electrode and headmodel units
% the output source model units should be consistent with the ones specified here (and not with elec or headmodel)
cfg.grid.resolution = 3;
cfg.grid.unit = 'cm';

elecunit      = {'mm', 'cm', 'm', 'inch'};
headmodelunit = {'mm', 'cm', 'm', 'feet'};

for i=1:numel(elecunit)
  for j=1:numel(headmodelunit)
    stopwatch = tic;
    
    fprintf('=================================================================================\n');
    fprintf('electrodes are expressed in %s, the headmodel is expressed in %s\n', elecunit{i}, headmodelunit{j});
    fprintf('=================================================================================\n');
    cfg.elec      = ft_convert_units(elec,      elecunit{i});
    cfg.headmodel = ft_convert_units(headmodel, headmodelunit{j});
    source{i,j} = ft_dipolefitting(cfg, timelock_data);
    
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
    if norm(source{i,j}.dip.pos - source{1,1}.dip.pos)>0.1 % cm
      error('difference in dipole fit results, please check the handling of geometrical units');
    end
  end
end

