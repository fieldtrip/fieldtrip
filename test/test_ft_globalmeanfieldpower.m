function test_ft_globalmeanfieldpower

% MEM 1500mb
% WALLTIME 00:10:00

% TEST test_ft_globalmeanfield
% TEST ft_singleplotER
load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2630/avg_tms_clean.mat'));

cfg = [];
cfg.channel = 'all';
cfg.method = 'power'; % amplitude/power

gmfp = ft_globalmeanfieldpower(cfg, avg_tms_clean);


cfg.method = 'amplitude';
gmfa = ft_globalmeanfieldpower(cfg, avg_tms_clean);

% plot
figure;
cfg = [];
ft_singleplotER(cfg, gmfa, gmfp);

end

