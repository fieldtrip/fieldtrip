function test_bug1502

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_checkconfig

% this test script should be explicit about ft_default
global ft_default
ft_default = [];

cfg   = [];
cfg.a = 1;
cfg.b = 2;
cfg.c = 3;

cfg = ft_checkconfig(cfg, 'allowed', {'a', 'b', 'c'});

% this field should not be allowed
cfg.d = 4;

try
  cfg = ft_checkconfig(cfg, 'allowed', {'a', 'b', 'c'});
  ok = false;
catch
  ok = true;
end

assert(ok, 'ft_checkconfig did not parse the allowed option correctly');
