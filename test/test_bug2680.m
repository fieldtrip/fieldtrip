function test_bug2680

% TEST ft_redefinetrial ft_timelockanalysis

% WALLTIME 00:10:00
% MEM 1gb

% create some test data, cut it into segments with ft_redefinetrial, and
% then call ft_timelockanalysis

data          = [];
data.label    = {'chan01';'chan02'};
data.trial{1} = randn(2,10000);
data.time{1}  = (0:9999)./1000;

cfg        = [];
cfg.length = 1;
datacut    = ft_redefinetrial(cfg, data);

cfg              = [];
try,
  tlck0            = ft_timelockanalysis(cfg, datacut);
  error('this should result in an error');
catch
  % this is OK
end
  
cfg.vartrllength = 1;
tlck1            = ft_timelockanalysis(cfg, datacut);

cfg.vartrllength = 2;
tlck2            = ft_timelockanalysis(cfg, datacut);
