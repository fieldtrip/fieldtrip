function test_bug2685

% MEM 4gb
% WALLTIME 00:10:00

% TEST test_bug2686
% TEST ft_scalpcurrentdensity ft_fetch_sens

%% load data
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2685/bug2685.mat'));

%% scalp current density
cfg                 = [];
cfg.method          = 'spline';
cfg.elec            = ERP_standard.elec;
%cfg.elec            = % Do we have electrode positions?

scd   = ft_scalpcurrentdensity(cfg, ERP_standard);
end
