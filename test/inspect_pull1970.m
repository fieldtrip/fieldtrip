function inspect_pull1970

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_databrowser

cfg = [];
cfg.dataset = dccnpath('/home/common/matlab/fieldtrip/data/ftp/test/ctf/Subject01.ds');
cfg = ft_definetrial(cfg);
cfg.trl = cfg.trl(1:20,:);
cfg.trl(:,3) = cfg.trl(:,1)-1;
cfg.demean = 'yes';
cfg.channel = 'MEG';
data = ft_preprocessing(cfg);

%%

cfg = [];
cfg.plotevents = 'no';
cfg.viewmode = 'vertical';
cfg.artfctdef.other.artifact = zeros(0,2);
cfg.artfctdef.blink.artifact = zeros(0,2);
cfg.artfctdef.movement.artifact = zeros(0,2);
ft_databrowser(cfg, data)

%%
% all viewmodes can work with component data

cfg = [];
cfg.method = 'pca';
comp = ft_componentanalysis(cfg, data);

%%

cfg = [];
cfg.plotevents = 'no';
cfg.viewmode = 'vertical';
ft_databrowser(cfg, comp)

%%

cfg = [];
cfg.plotevents = 'no';
cfg.viewmode = 'butterfly';
ft_databrowser(cfg, comp)

%%

cfg = [];
cfg.plotevents = 'no';
cfg.viewmode = 'component';
ft_databrowser(cfg, comp)
