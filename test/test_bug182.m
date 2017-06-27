function test_bug182

% MEM 1500mb
% WALLTIME 00:10:00

% TEST ft_preprocessing ft_componentanalysis ft_rejectcomponent ft_componentbrowser ft_databrowser

% this script addresses bug 182.
% applying ft_componentanalysis and reconstructing the data with
% ft_rejectcomponent should have the possibility to contain a
% grad-structure which is balanced according to the mixing matrix

cd(dccnpath('/home/common/matlab/fieldtrip/data/test/latest/raw/meg/'));
load('preproc_ctf151.mat');

% the demean is essential to have an equal output for datanew1 and datanew2
cfg = [];
cfg.channel = {'MEG'};
cfg.demean  = 'yes';
data = ft_preprocessing(cfg, data);

cfg = [];
cfg.method = 'runica';
cfg.channel = {'MEG'};
cfg.runica.maxsteps = 100;
comp = ft_componentanalysis(cfg, data);

% visualize
cfg = [];
cfg.component = 1:20;
cfg.layout = 'CTF151.lay';
cfg.viewmode = 'component';
ft_databrowser(cfg, comp);

% now we can call ft_rejectcomponent in two ways,
% -reconstruct the data using comp only
cfg = [];
cfg.component = 7;
datanew1 = ft_rejectcomponent(cfg, comp);

% -reconstruct the data using the mixing matrix from comp and the original
% data
datanew2 = ft_rejectcomponent(cfg, comp, data);

figure;plot(datanew1.trial{1}(1,:));hold on;plot(datanew2.trial{1}(1,:),'r');

% this should reconstruct the original data, including the original grad
cfg = [];
cfg.component = [];
datanew3 = ft_rejectcomponent(cfg, comp);
