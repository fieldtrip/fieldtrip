function test_bug3293

% WALLTIME 00:10:00
% MEM 1gb

%%

elec = [];
elec.label = {'1', '2', '3', '4'};
elec.elecpos = randn(4,3);
elec.chanpos = elec.elecpos;
elec.tra = eye(4);

try
  cfg = [];
  ft_selectdata(cfg, elec);
  failed = false;
catch
  failed = true;
end
assert(failed, 'it should have failed');


