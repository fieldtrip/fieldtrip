function test_bug2059

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_topoplotER ft_channelselection

load(dccnpath('/home/common/matlab/fieldtrip/data/test/bug2059.mat'));

close all

cfg = [];
cfg.channel = {'4','5','6','8','9','10','11','13','14','15','19','20','21','24','25','26'};
cfg.layout = lay_baby;
cfg.xlim = [0.7 1.2];
cfg.zlim = [-6 6];
% ft_topoplotER(cfg, GA_loc_diff);

cfg.channel = 'all';
figure(1)
ft_topoplotER(cfg, GA_loc_diff);

cfg.channel = {'4','11','16','19','20','21'};
figure(2)
ft_topoplotER(cfg, GA_loc_diff);

% these figures should be different
assert(~isequalfigure(1, 2));

% the problem seemed to be in a buggy interaction between interactive and channel selection

cfg.interactive = 'no';
cfg.channel = {'4','11','16','19','20','21'};
figure(3)
ft_topoplotER(cfg, GA_loc_diff);


cfg.interactive = 'yes';
cfg.channel = {'4','11','16','19','20','21'};
figure(4)
ft_topoplotER(cfg, GA_loc_diff);

% these figures should be the same
assert(isequalfigure(3, 4));

%% specifying layout as a string

cfg.layout = 'easycapM14.lay';
cfg.interactive = 'no';
cfg.channel = {'4','11','16','19','20','21'};
figure(5)
ft_topoplotER(cfg, GA_loc_diff);


cfg.interactive = 'yes';
cfg.channel = {'4','11','16','19','20','21'};
figure(6)
ft_topoplotER(cfg, GA_loc_diff);

% read in layout directly
cfg.layout = ft_prepare_layout(cfg);
cfg.interactive = 'no';
cfg.channel = {'4','11','16','19','20','21'};
figure(7)
ft_topoplotER(cfg, GA_loc_diff);


cfg.interactive = 'yes';
cfg.channel = {'4','11','16','19','20','21'};
figure(8)
ft_topoplotER(cfg, GA_loc_diff);



% all these figures should be the same
assert(isequalfigure(5, 6));
assert(isequalfigure(5, 7));
assert(isequalfigure(5, 8));
assert(isequalfigure(6, 7));
assert(isequalfigure(6, 8));
assert(isequalfigure(7, 8));

