function test_bug2685

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_scalpcurrentdensity ft_fetch_sens
% DATA private

%% load data
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2685/bug2685.mat'));

%% scalp current density
cfg                 = [];
cfg.method          = 'spline';
cfg.elec            = ERP_standard.elec;
%cfg.elec            = % Do we have electrode positions?

scd   = ft_scalpcurrentdensity(cfg, ERP_standard);
end
