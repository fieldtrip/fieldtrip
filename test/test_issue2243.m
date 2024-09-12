function test_issue2243

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_redefinetrial
% DATA no

% Create dummy data
data = [];
data.trial = mat2cell(randn(3, 4004), 3, [1001 1001 1001 1001]);
data.time  = repmat({linspace(0, 1, 1001)}, 1, 4);
data.label = {'A', 'B', 'C'};

% This incorrectly shifts the second trial by 0.3 seconds
cfg        = [];
cfg.trials = [false, true, false, false];
cfg.offset = [-500, -300];
try
  data_out = ft_redefinetrial(cfg, data);
catch me
  assert(strcmp(me.message, 'inconsistent number of trials and offsets'));
end

% This throws an error because it looks for more offsets than it got
cfg        = [];
cfg.trials = [false, true, true, false];
cfg.offset = [-500, -300];
data_out = ft_redefinetrial(cfg, data);

% This throws an error because not enough toilims are provided
cfg        = [];
cfg.trials = [false, true, true, true];
cfg.toilim = [0, 0.3; 0.3, 0.5];
try
  data_out = ft_redefinetrial(cfg, data);
catch me
  assert(strcmp(me.message, 'inconsistent number of trials and toilims'));
end
