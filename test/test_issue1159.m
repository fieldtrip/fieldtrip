function test_issue1159

% WALLTIME 00:10:00
% MEM 2gb
% DEPENDENCY ft_senstype undobalancing

%%
% see https://github.com/fieldtrip/fieldtrip/issues/1159

load issue1159.mat

%%
% this works

ft_senstype(grad)

%%
% this originally fails

ft_senstype(rmfield(grad, 'type'))

%%

cfg = [];
cfg.dataset = '/home/common/matlab/fieldtrip/data/test/original/meg/neuromag306/raw.fif';
data = ft_preprocessing(cfg);

cfg = [];
cfg.numcomponent = 12;
cfg.method = 'pca';
cfg.channel = 'MEG';
comp = ft_componentanalysis(cfg, data);


ft_senstype(data.grad) % this works
ft_senstype(comp.grad) % this originally fails
