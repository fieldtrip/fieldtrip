function test_pull1331

% MEM 2gb
% WALLTIME 00:10:00
% DEPENDENCY ft_prepare_layout ft_preprocessing ft_channelrepair

load(fullfile(fileparts(which('ft_defaults')),'template/neighbours/ctf151_neighb.mat'));
lay = ft_prepare_layout(struct('layout',fullfile(fileparts(which('ft_defaults')),'template/layout/CTF151.lay')));

ds = 'Subject01.ds';
cfg = [];
cfg.dataset = ds;
data = ft_preprocessing(cfg);


cfg = [];
cfg.badchannel =     {'MLC11','MLC12','MLC13' };

cfg.neighbours = neighbours;
cfg.layout = lay;
data_clean = ft_channelrepair(cfg,data);

