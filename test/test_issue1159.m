function test_issue1159

% WALLTIME 00:10:00
% MEM 1gb
% DEPENDENCY ft_senstype undobalancing
% DATA private

%%
% see https://github.com/fieldtrip/fieldtrip/issues/1159

load(dccnpath('/project/3031000.02/test/issue1159.mat'))

%%
% this works

ft_senstype(grad)

%%
% this originally fails

ft_senstype(rmfield(grad, 'type'))

%%
% the following is a generalization and also applies to https://github.com/fieldtrip/fieldtrip/pull/980

cfg = [];
cfg.dataset = dccnpath('/project/3031000.02/test/original/meg/neuromag306/raw.fif');
cfg.trl = [1 7500 0]; % 30 seconds @250Hz
cfg.continuous = 'yes';
data = ft_preprocessing(cfg);

% hmm, there is another (not related) error, see https://github.com/fieldtrip/fieldtrip/issues/1162
% it relates to channels that are all-zero being tossed out of the montage
% this is a quick fix
data.trial{1}(1:3,:) = randn(size(data.trial{1}(1:3,:)));

cfg = [];
cfg.numcomponent = 12;
cfg.method = 'pca'; % the type of unmixing does not matter, this is the fastest
cfg.channel = {'MEG'};
comp = ft_componentanalysis(cfg, data);

% this originally failed
cfg = [];
cfg.component = [1 2 3];
clean = ft_rejectcomponent(cfg, comp, data);

ft_senstype(data.grad) % this works
ft_senstype(comp.grad) % this originally fails

%%
% try the same thing with a CTF dataset

cfg = [];
cfg.dataset = dccnpath('/project/3031000.02/external/download/test/ctf/Subject01.ds');
cfg.trl = [1 9000 0]; % 30 seconds @300Hz
cfg.continuous = 'yes';
data = ft_preprocessing(cfg);

cfg = [];
cfg.numcomponent = 12;
cfg.method = 'pca'; % the type of unmixing does not matter, this is the fastest
cfg.channel = 'MEG';
comp = ft_componentanalysis(cfg, data);

% this originally failed
cfg = [];
cfg.component = [1 2 3];
clean = ft_rejectcomponent(cfg, comp, data);

ft_senstype(data.grad) % this works
ft_senstype(comp.grad) % this originally fails
