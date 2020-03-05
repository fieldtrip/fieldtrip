function inspect_issue1216

% WALLTIME 00:10:00
% MEM 3gb
% DEPENDENCY ft_hastoolbox ft_postamble_hastoolbox

% this should be executed on MATLAB version R2016a, since it deals with the compat/matlablt2016b directory

load(dccnpath('/home/common/matlab/fieldtrip/data/test/issue1216.mat'));

%% Run ICA
% Perform independent component analysis
cfg = [];
cfg.method = 'runica'; % this is the default and uses the implementation from EEGLAB
comp = ft_componentanalysis(cfg, data_NoBigArt);

% for an unclear reason the code above returns a complex comp.topo and comp.trial on R2016a
% the issue does not relate to this, hence I am recomputing with PCA

cfg = [];
cfg.method = 'pca'; % this is the default and uses the implementation from EEGLAB
comp = ft_componentanalysis(cfg, data_NoBigArt);


%% Identify artefacts
% Visualize the time course of the artefacts
cfg = [];
cfg.viewmode = 'component';
cfg.continuous = 'no'; % view the data in it's segments
cfg.channel = 'all'; % view all the component in one screen
cfg.layout = lay; % layout for DCC acticap
ft_databrowser(cfg,comp);

%% Visualize the components
% This is were the error occurs, the error states that there is an
% undefined function 'newline' which is located in the Matlablt2016b
% folder. However, this folder is removed from the path before executing
% the function
figure
cfg = [];
cfg.component = 1:20;
cfg.layout = lay'; % specify the layout file that should be used for plotting
cfg.comment = 'no';
ft_topoplotIC(cfg, comp);

