function test_bug3433

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_heartrate ft_respiration

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/bug3433'));

%%
cfg =[];
cfg.dataset = '006_3013065.02_rest2.vhdr';
data = ft_preprocessing(cfg);

%%

cfg = [];
cfg.toilim = [0 60];
data = ft_redefinetrial(cfg, data);

%%

cfg = [];
cfg.resamplefs = 500;
data = ft_resampledata(cfg, data);

%%

cfg = [];
cfg.channel = 'HR';
cfg.threshold = 0.7;
heartrate = ft_heartrate(cfg, data);


figure
plot(data.time{1}, ft_preproc_standardize(data.trial{1}(1,:)));
hold on
plot(heartrate.time{1}, heartrate.trial{1}(3,:), 'r');
plot(heartrate.time{1}, heartrate.trial{1}(2,:), 'g');

%%

cfg = [];
cfg.channel = 'Resp';
respiration = ft_respiration(cfg, data);

