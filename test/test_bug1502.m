function test_bug1502

% MEM 1gb
% WALLTIME 00:10:00
% DEPENDENCY ft_checkconfig
% DATA no

cfg   = [];
cfg.a = 1;
cfg.b = 2;
cfg.c = 3;
cfg.checkconfig = 'pedantic'; % this field pertains to error handling and is not checked itself

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
