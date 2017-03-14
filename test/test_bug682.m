function test_bug682

% MEM 1500mb
% WALLTIME 00:10:00

% TEST channelposition ft_datatype_sens yokogawa2grad ft_prepare_layout

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/original/meg/yokogawa64'));

cfg             = [];
cfg.dataset     = '2011_01_28_0354_ME053_AEF.con'; % to be defined beforehand
cfg.demean      = 'no';
cfg.channel     = 'AG*';
% magnetometers: cfg.channel     = 'M*';
% planar gradiometers: cfg.channel     = 'PG*';
% all: cfg.channel     = 'all';

data            = ft_preprocessing(cfg);
% data.grad.label contains only valid channels, bad channels marked
% by the user in the Yokogawa MegLaboratory software are excluded
% (MegLaboratoy: Open original file, then menu: Edit ->  Channel Info,
% select channel, change Channel Type to Nullchannel; Then save file
% using a new name)

% plot sensor geometry
figure('name','sensor');
nr_chan          = size(data.grad.chanpos,1);
plot3(data.grad.chanpos(1:nr_chan,1), data.grad.chanpos(1:nr_chan,2),data.grad.chanpos(1:nr_chan,3),'d')
xlim([-15,15]); zlim([-15,15]); ylim([-15,15]);

hold on
ft_plot_sens(data.grad);

% plot sensor layout
cfg         = [];
cfg.layout  = ft_prepare_layout([], data);
figure('name','layout');
ft_layoutplot(cfg, data);
